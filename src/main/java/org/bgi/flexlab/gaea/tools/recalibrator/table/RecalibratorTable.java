package org.bgi.flexlab.gaea.tools.recalibrator.table;

import java.util.ArrayList;
import java.util.List;

import org.bgi.flexlab.gaea.tools.mapreduce.realigner.RecalibratorOptions;
import org.bgi.flexlab.gaea.tools.recalibrator.RecalibratorDatum;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.Covariate;
import org.bgi.flexlab.gaea.tools.recalibrator.covariate.CovariateUtil;
import org.bgi.flexlab.gaea.util.EventType;
import org.bgi.flexlab.gaea.util.NestedObjectArray;
import org.bgi.flexlab.gaea.util.NestedObjectArray.Leave;

import htsjdk.samtools.SAMFileHeader;

public class RecalibratorTable {
	public enum Type {
		READ_GROUP_TABLE(0), QUALITY_SCORE_TABLE(1), OPTIONAL_COVARIATE_TABLES_START(2);

		public final int index;

		private Type(final int index) {
			this.index = index;
		}
	}

	@SuppressWarnings("rawtypes")
	private NestedObjectArray[] tables = null;

	public RecalibratorTable(final Covariate[] covariates, int readGroupNumber) {
		tables = new NestedObjectArray[covariates.length];

		int maxQualityScore = covariates[Type.QUALITY_SCORE_TABLE.index].maximumKeyValue() + 1;
		int eventSize = EventType.values().length;

		tables[Type.READ_GROUP_TABLE.index] = new NestedObjectArray<RecalibratorDatum>(readGroupNumber, eventSize);
		tables[Type.QUALITY_SCORE_TABLE.index] = new NestedObjectArray<RecalibratorDatum>(readGroupNumber,
				maxQualityScore, eventSize);
		for (int i = Type.OPTIONAL_COVARIATE_TABLES_START.index; i < covariates.length; i++)
			tables[i] = new NestedObjectArray<RecalibratorDatum>(readGroupNumber, maxQualityScore,
					covariates[i].maximumKeyValue() + 1, eventSize);
	}

	public static RecalibratorTable build(RecalibratorOptions option, SAMFileHeader header) {
		RecalibratorTable recalibratorTables = new RecalibratorTable(CovariateUtil.initializeCovariates(option, header),
				header.getReadGroups().size());

		return recalibratorTables;
	}

	@SuppressWarnings("unchecked")
	public NestedObjectArray<RecalibratorDatum> getTable(int index) {
		return (NestedObjectArray<RecalibratorDatum>) tables[index];
	}

	public NestedObjectArray<RecalibratorDatum> getTable(Type type) {
		return getTable(type.index);
	}

	public int length() {
		return tables.length;
	}

	public ArrayList<String> valueStrings() {
		ArrayList<String> arrays = new ArrayList<String>();
		for (int i = 0; i < tables.length; i++) {
			@SuppressWarnings("unchecked")
			NestedObjectArray<RecalibratorDatum> table = tables[i];
			List<Leave> leaves = table.getAllLeaves();
			for (Leave leave : leaves) {
				arrays.add(leave.toString(i));
			}
		}
		return arrays;
	}
}