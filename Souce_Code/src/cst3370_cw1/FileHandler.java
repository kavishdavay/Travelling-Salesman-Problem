package cst3370_cw1;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

/**
 * File handler to read files and return integer values in a 2D array.
 * Has validation to catch incorrect file format or missing file.
 * 
 * @version 1.0 06 Dec 2023
 * @author M00789089
 * 
 */
public class FileHandler {

	/**
	 * Method to read a file and output a 2D array with file contents.
	 * 
	 * @param filePath The file path
	 * @return 2D array of integers
	 */
	public static int[][] read(String filePath) {
		
		/* Try to open file using file path */
		try {
			
			/* Initialise scanner with file path */
			Scanner reader = new Scanner(new File(filePath));
			
			/* Initialise row and column count */
			int rowNum = 0;
			int columnNum = 0;
			
			/* Execute for each line in file */
			while (reader.hasNextLine()) {
				
				/* Increment row count */
				rowNum++;
				
				/* 
				 * Determine number of columns in file.
				 * Assumes that each column is separated by a tab.
				 * Trim is used to account for spaces before and after values in a column.
				 */
				String[] data = reader.nextLine().trim().split("\\s+");
				columnNum = data.length;
			}
			
			/* Checks if file has correct number of columns */
			if(columnNum == 3) {
				
				/* Re-initialise scanner */
				reader = new Scanner(new File(filePath));
				
				/* Initialise integer 2D array with row and column count */
				int[][] integerArray = new int[rowNum][columnNum];
				
				/* Go through each row and column */
				for (int i = 0; i < rowNum; i++) {
					for (int j = 0; j < columnNum; j++) {
						
						/* Check for integer at index position */
						if (reader.hasNextInt()) {
							
							/* Add integer to 2D array */
							integerArray[i][j] = reader.nextInt();
							
						} else {
							
							/* Throw an error if a non-integer is found */
							throw new IllegalStateException("Non-integer value found");
						}
					}
				}
				
				/* Close the scanner */
				reader.close();
				
				return integerArray;
				
			} else {
				
				/* Throw an error if file doesn't have 3 columns */
				throw new IllegalStateException("Incorrect file format");
				
			}

			
		} catch (FileNotFoundException e) {
			
			/* Throw error if file not found */
            throw new RuntimeException("File not found: " + filePath);
		}
    }
	
}
