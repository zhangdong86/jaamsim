/*
 * JaamSim Discrete Event Simulation
 * Copyright (C) 2005-2011 Ausenco Engineering Canada Inc.
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
 */
package com.sandwell.JavaSimulation;

import java.util.ArrayList;


	/**
	 *
	 *
	 */
public class Util {

	/**
	 * Expects a StringVector of one of two forms:
	 * 1.  { entry entry entry } { entry entry entry }
	 * 2.  entry entry entry
	 * If format 1, returns a vector of stringvectors without braces.
	 * if format 2, returns a vector of a stringvector, size1.
	 * @param data
	 * @return
	 */
	public static ArrayList<StringVector> splitStringVectorByBraces(StringVector data) {
		ArrayList<StringVector> newData = new ArrayList<StringVector>();
		for (int i=0; i < data.size(); i++) {

			//skip over opening brace if present
			if (data.get(i).equals("{") )
				continue;

			StringVector cmd = new StringVector();

			//iterate until closing brace, or end of entry
			for (int j = i; j < data.size(); j++, i++){
				if (data.get(j).equals("}"))
					break;

				cmd.add(data.get(j));
			}

			//add to vector
			newData.add(cmd);
		}

		return newData;
	}
}
