/*******************************************************************************
 * Copyright (c) 2017, BGI-Shenzhen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 *******************************************************************************/
package org.bgi.flexlab.gaea.data.mapreduce.partitioner;

import org.apache.hadoop.mapreduce.Partitioner;
import org.bgi.flexlab.gaea.data.mapreduce.writable.WindowsBasedWritable;

public class WindowsBasedPartitioner<T> extends Partitioner<WindowsBasedWritable, T> {

	@Override
	public int getPartition(WindowsBasedWritable key, T v, int numPartitioner) {
		int hashcode = key.partition() ;
		return Math.abs(hashcode) % numPartitioner;
	}
}
