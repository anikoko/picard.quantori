/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.AbstractLocusInfo;
import htsjdk.samtools.util.AbstractLocusIterator;
import htsjdk.samtools.util.AbstractRecordAndOffset;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.Log;
import picard.filter.CountingFilter;
import picard.filter.CountingPairedFilter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Semaphore;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.LongStream;

/**
 * Implementation of {@link picard.analysis.WgsMetricsProcessor} that gets input data from a given iterator
 * and processes it with a help of collector
 * @author Mariia_Zueva@epam.com, EPAM Systems, Inc. <www.epam.com>
 */
public class WgsMetricsProcessorImpl<T extends AbstractRecordAndOffset> implements WgsMetricsProcessor {
    private static final int BATCH_SIZE = 1000;
    /**
     * Source of input data
     */
    private final AbstractLocusIterator<T, AbstractLocusInfo<T>> iterator;
    /**
     * Accumulates the data from iterator
     */
    private final AbstractWgsMetricsCollector<T> collector;
    /**
     * ReferenceWalker for a processed reference sequence
     */
    private final ReferenceSequenceFileWalker refWalker;
    /**
     * Logger for the progress of work
     */
    private final ProgressLogger progress;

    private final Log log = Log.getInstance(WgsMetricsProcessorImpl.class);

    /**
     * @param iterator  input {@link htsjdk.samtools.util.AbstractLocusIterator}
     * @param refWalker over processed reference file
     * @param collector input {@link picard.analysis.AbstractWgsMetricsCollector}
     * @param progress  logger
     */
    public WgsMetricsProcessorImpl(AbstractLocusIterator<T, AbstractLocusInfo<T>> iterator,
                                   ReferenceSequenceFileWalker refWalker,
                                   AbstractWgsMetricsCollector<T> collector,
                                   ProgressLogger progress) {
        this.iterator = iterator;
        this.collector = collector;
        this.refWalker = refWalker;
        this.progress = progress;
    }

    /**
     * Method gets the data from iterator for each locus and processes it with the help of collector.
     */
    @Override
    public void processFile() {
        AtomicLong counter = new AtomicLong(0L);

        List<AbstractLocusInfo<T>> batch = new ArrayList<>(BATCH_SIZE);

        ExecutorService service = Executors.newSingleThreadExecutor();

        Semaphore sem = new Semaphore(1);

        while (iterator.hasNext()) {
            final AbstractLocusInfo<T> info = iterator.next();
            batch.add(info);

            if (batch.size() == BATCH_SIZE) {
                submitInfos(batch, service, counter, sem);
                batch = new ArrayList<>(BATCH_SIZE);
            }

        }

        if (!batch.isEmpty()){
            submitInfos(batch, service, counter, sem);
            batch = new ArrayList<>(BATCH_SIZE);
        }
        // check that we added the same number of bases to the raw coverage histogram and the base quality histograms
        final long sumBaseQ = Arrays.stream(collector.unfilteredBaseQHistogramArray).sum();
        final long sumDepthHisto = LongStream.rangeClosed(0, collector.coverageCap).map(i -> (i * collector.unfilteredDepthHistogramArray[(int) i])).sum();
        if (sumBaseQ != sumDepthHisto) {
            log.error("Coverage and baseQ distributions contain different amount of bases!");
        }
    }

    private void submitInfos(List<AbstractLocusInfo<T>> batch, ExecutorService service, AtomicLong counter, Semaphore sem) {
        List<AbstractLocusInfo<T>> infos = batch;
        service.submit(() -> {
            for (AbstractLocusInfo<T> absInfo : infos) {
                final ReferenceSequence ref = refWalker.get(absInfo.getSequenceIndex());
                boolean referenceBaseN = collector.isReferenceBaseN(absInfo.getPosition(), ref);
                collector.addInfo(absInfo, ref, referenceBaseN);
                if (referenceBaseN) {
                    continue;
                }

                progress.record(absInfo.getSequenceName(), absInfo.getPosition());
                if (collector.isTimeToStop(counter.incrementAndGet())) {
                    break;
                }
                collector.setCounter(counter.get());

            }
            sem.release();
        });
    }

    @Override
    public void addToMetricsFile(MetricsFile<WgsMetrics, Integer> file,
                                 boolean includeBQHistogram,
                                 CountingFilter dupeFilter,
                                 CountingFilter adapterFilter,
                                 CountingFilter mapqFilter,
                                 CountingPairedFilter pairFilter) {
        collector.addToMetricsFile(file, includeBQHistogram, dupeFilter, adapterFilter, mapqFilter, pairFilter);
    }
}
