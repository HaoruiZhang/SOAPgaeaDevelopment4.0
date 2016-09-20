package org.bgi.flexlab.gaea.data.structure.pileup;

import java.util.List;
import java.util.Map;

import org.bgi.flexlab.gaea.data.structure.bam.GaeaSamRecord;
import org.bgi.flexlab.gaea.data.structure.location.GenomeLocation;

public class PileupImpl extends AbstractPileup<PileupImpl, PileupElement> {
	public PileupImpl(GenomeLocation loc) {
		super(loc);
	}

	public PileupImpl(GenomeLocation loc, List<GaeaSamRecord> reads,
			List<Integer> offsets) {
		super(loc, reads, offsets);
	}

	public PileupImpl(GenomeLocation loc, List<GaeaSamRecord> reads, int offset) {
		super(loc, reads, offset);
	}

	public PileupImpl(GenomeLocation loc, List<PileupElement> pileupElements) {
		super(loc, pileupElements);
	}

	public PileupImpl(GenomeLocation loc,
			Map<String, PileupImpl> pileupElementsBySample) {
		super(loc, pileupElementsBySample);
	}

	/**
	 * Optimization of above constructor where all of the cached data is
	 * provided
	 */
	public PileupImpl(GenomeLocation loc, List<PileupElement> pileup, int size,
			int nDeletions, int nMQ0Reads) {
		super(loc, pileup, size, nDeletions, nMQ0Reads);
	}

	protected PileupImpl(GenomeLocation loc,
			PileupElementTracker<PileupElement> tracker) {
		super(loc, tracker);
	}

	@Override
	protected PileupImpl createNewPileup(GenomeLocation loc,
			PileupElementTracker<PileupElement> tracker) {
		return new PileupImpl(loc, tracker);
	}

	@Override
	protected PileupElement createNewPileupElement(final GaeaSamRecord read,
			final int offset, final boolean isDeletion,
			final boolean isBeforeDeletion, final boolean isAfterDeletion,
			final boolean isBeforeInsertion, final boolean isAfterInsertion,
			final boolean isNextToSoftClip) {
		return new PileupElement(read, offset, isDeletion, isBeforeDeletion,
				isAfterDeletion, isBeforeInsertion, isAfterInsertion,
				isNextToSoftClip, null, 0);
	}

	@Override
	protected PileupElement createNewPileupElement(final GaeaSamRecord read,
			final int offset, final boolean isDeletion,
			final boolean isBeforeDeletion, final boolean isAfterDeletion,
			final boolean isBeforeInsertion, final boolean isAfterInsertion,
			final boolean isNextToSoftClip, final String nextEventBases,
			final int nextEventLength) {
		return new PileupElement(read, offset, isDeletion, isBeforeDeletion,
				isAfterDeletion, isBeforeInsertion, isAfterInsertion,
				isNextToSoftClip, nextEventBases, nextEventLength);
	}
}
