package cst3370_cw1;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Used to store genes and has a fitness score.
 * Includes a method to randomly shuffle genes for the initial population.
 * 
 * @version 1.0 06 Dec 2023
 * @author M00789089
 * 
 */
public class Chromosome {
	private int[] genes;						// Stores city IDs
	private double fitnessScore = -1;			// Initialise fitness score
	
	/**
	 * Initialise chromosome with integer array of genes.
	 * 
	 * @param genes
	 */
	public Chromosome(int[] genes) {
		this.genes = genes;
	}
	
	/**
	 * Get genes from chromosome.
	 * 
	 * @return Integer array of genes
	 */
	public int[] getChromosome() {
		return this.genes;
	}
	
	/**
	 * Get fitness score
	 * 
	 * @return fitness score
	 */
	public double getFitnessScore() {
		return this.fitnessScore;
	}
	
	/**
	 * Get a gene at a specific position in the chromosome.
	 * 
	 * @param genePosition
	 * @return Gene
	 */
	public int getGene(int genePosition) {
		return this.genes[genePosition];
	}
	
	/**
	 * Set a gene at a specific position in the chromosome.
	 * 
	 * @param index Position in chromosome
	 * @param gene
	 */
	public void setGene(int index, int gene) {
		this.genes[index] = gene;
	}
	
	/**
	 * Set fitness score
	 * 
	 * @param fitnessScore
	 */
	public void setFitnessScore (double fitnessScore) {
		this.fitnessScore = fitnessScore;
	}
	
	/**
	 * Method to randomly shuffle genes in the chromosome.
	 */
	public void shuffleGenes() {
		
		/* Creating ArrayList to store shuffled genes */
		List<Integer> genesList = new ArrayList<>();
		
		/* Go through each gene in the chromosome and add it to the ArrayList */
		for(int gene : this.genes) {
			genesList.add(gene);
		}
		
		/* Use of Collections to randomly shuffle genes in the ArrayList */
		Collections.shuffle(genesList);
		
		/* Replacing original genes with shuffled genes */
		for(int i = 0; i < this.genes.length; i++) {
			this.genes[i] = genesList.get(i);
		}
	}
	
}
