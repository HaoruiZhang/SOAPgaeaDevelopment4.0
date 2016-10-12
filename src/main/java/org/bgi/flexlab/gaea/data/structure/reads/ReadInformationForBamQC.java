/**
 * Copyright (c) 2011, BGI and/or its affiliates. All rights reserved.
 * BGI PROPRIETARY/CONFIDENTIAL. Use is subject to license terms.
 */
package org.bgi.flexlab.gaea.data.structure.reads;

import org.bgi.flexlab.gaea.util.CigarState;
import org.bgi.flexlab.gaea.util.ParseSAMBasic;
import org.bgi.flexlab.gaea.util.ParseSAMInterface;


/**
 * 序列比对结果类
 * @author ZhangYong
 *
 */
public class ReadInformationForBamQC extends SAMInformationBasic {
	/**
	 * Read ID
	 */
	private String id;

	/**
	 * Read比对到参考基因组上的次数
	 */
	private int bestHitCount;
	
	/**
	 * insert Size
	 */
	private int insertSize;
	
	/**
	 * 
	 */
	private boolean isrepeat = false;
	
	/**
	 * read end on ref, reads seq length in sam, reads base count on ref
	 */
	private int[] lenValue;
	
	/**
	 * rg index
	 */
	private int rgIndex;
	
	@Override
	protected void parseOtherInfo(String[] alignmentArray) {

		bestHitCount = ParseSAMBasic.parseBestHitCount(alignmentArray);
		
		insertSize = ParseSAMBasic.parseInsertSize(alignmentArray);
		
		lenValue = ParseSAMBasic.parseCigar(position, cigarState);
	}
	
	public boolean parseBAMQC(String value) {
		// 该行为空时，不进行解析
		if (value.isEmpty()) {
			return false;
		}

		// 以制表符分割SAM比对文件的每行字符串
		String[] alignmentArray = value.split("\t");
		
		flag = Integer.parseInt(alignmentArray[0]);
		
		readSequence = alignmentArray[1];
		
		insertSize = Integer.parseInt(alignmentArray[2]);
		
		position = Integer.parseInt(alignmentArray[3]);
		
		if(position < 0 ) {
			return false;
		}
		
		cigarString = alignmentArray[4];
		cigarState = new CigarState();
		cigarState.parseCigar(cigarString);
		
		bestHitCount = Integer.parseInt(alignmentArray[5]);
		
		if(Integer.parseInt(alignmentArray[6]) == 1) {
			isrepeat = true;
		}
		
		mappingQual = Short.parseShort(alignmentArray[7]);
		if(alignmentArray.length >= 9) {
			rgIndex = Integer.parseInt(alignmentArray[8]);
		}
		if(alignmentArray.length >= 10) {
			qualityString = alignmentArray[9];
		}
		lenValue = ParseSAMBasic.parseCigar(position, cigarState);
		return true;
	}

	/**
	 * 获取Read比对到参考基因组上的最佳比对次数
	 * @return int
	 */
	public int getBestHitCount() {
		return bestHitCount;
	}

	/**
	 * 判断Read比对到参考基因组上的次数是否为1
	 * @return true or false
	 */
	public boolean isUniqueAlignment() {
		return (bestHitCount == 1);
	}

	/**
	 * 获取Read ID
	 * @return
	 */
	public String getId() {
		return id;
	}

	/**
	 * @return the end
	 */
	public int getEnd() {
		return lenValue[0];
	}

	public int getLength() {
		return lenValue[1];
	}

	/**
	 * @return the insertSize
	 */
	public int getInsertSize() {
		return insertSize;
	}

	/**
	 * @return the isrepeat
	 */
	public boolean isIsrepeat() {
		return isrepeat;
	}

	/**
	 * @param isrepeat the isrepeat to set
	 */
	public void setIsrepeat(boolean isrepeat) {
		this.isrepeat = isrepeat;
	}

	/**
	 * @return the baseCount
	 */
	public int getBaseCount() {
		return lenValue[2];
	}

	/**
	 * @return the rgIndex
	 */
	public int getRgIndex() {
		return rgIndex;
	}

	/**
	 * @param rgIndex the rgIndex to set
	 */
	public void setRgIndex(int rgIndex) {
		this.rgIndex = rgIndex;
	}
}
