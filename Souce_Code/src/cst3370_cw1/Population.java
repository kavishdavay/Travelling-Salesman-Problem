package cst3370_cw1;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * This class stores a population of chromosomes.
 * The chromosomes are used to calculate the average fitness score of the population.
 * 
 * @version 1.0 06 Dec 2023
 * @author M00789089
 * 
 */
public class Population {

	private ArrayList<Chromosome> population;	// ArrayList of chromosome objects.
	private double averageFitness = -1;			// Initialising average fitness.
	
	/**
	 * Initialise an empty population with a specified size.
	 * 
	 * @param populationSize The size of the population
	 */
	public Population(int populationSize) {
		this.population = new ArrayList<>(populationSize);
	}
	
	/**
	 * Method to initialise first population.
	 *  
	 * @param populationSize Size of population
	 * @param cities Array containing city objects
	 */
	public Population(int populationSize, City cities[]) {
		this.population = new ArrayList<>(populationSize);		// Create new population
		int[] cityIDs = new int[cities.length];				// Array to store city IDs
		
		/* 
		 * Go through each city object in the cities array.
		 * Get the city ID and add it to the cityIDs array.
		 */
		for (int i = 0; i < cities.length; i++) {
			cityIDs[i] = cities[i].getCityID();
		}
		
		/* Go through the population to add chromosomes */
		for(int chromosomePos = 0; chromosomePos < populationSize; chromosomePos++) {
			
			/* Create a new chromosome with the cityIDs array */
			Chromosome chromosome = new Chromosome(cityIDs);
			
			/* Calling the shuffle method on the chromosome*/
			chromosome.shuffleGenes();
			
			/* Add the shuffled chromosome to the population */
			this.population.add(chromosomePos, chromosome);
		}
	}
	
	/**
	 * Get the population of chromosomes.
	 * 
	 * @return ArrayList of chromosomes
	 */
	public ArrayList<Chromosome> getPopulation() {
		return this.population;
	}
	
	/**
	 * Get a chromosome at a specific position in the population.
	 * 
	 * @param index Position of chromosome
	 * @return Chromosome object
	 */
	public Chromosome getChromosome(int index) {
		return this.population.get(index);
	}
	
	/**
	 * Get the average fitness of the population
	 * 
	 * @return double averageFitness
	 */
	public double getAverageFitness() {
		
		/* Check if average fitness has been calculated */
		if(this.averageFitness > -1) {
			return this.averageFitness;
		} else {
			setAverageFitness();			// Calling method to calculate average fitness
			return this.averageFitness;
		}
	}
	
	/**
	 * Calculate and set the average fitness of the population.
	 */
	public void setAverageFitness() {
		double fitnessSum = 0;				// Initialising total fitness
		
		/* 
		 * Go through each chromosome in the population,
		 * get its fitness score and add it to fitnessSum.
		 */
		for(Chromosome chromosome : this.population) {
			fitnessSum += chromosome.getFitnessScore();
		}
		
		/* Calculate average fitness by dividing fitness sum by population size*/
		double averageFitness = fitnessSum / this.population.size();
		
		/* Set population object average fitness */
		this.averageFitness = averageFitness;
	}
	
	/**
	 * Get chromosome at a position after sorting population.
	 * 
	 * @param index Chromosome index
	 * @return Chromosome object
	 */
	public Chromosome getBestChromosome(int index) {
		
		/* Call method to sort population by fitness score in descending order*/
		sortPopulation(this.population);
		
		/* Return chromosome at specified index */
		return this.population.get(index);
	}
	
	/**
	 * Sort the population of chromosomes by fitness score in descending order.
	 * 
	 * @param population ArrayList of chromosomes
	 */
	public static void sortPopulation(ArrayList<Chromosome> population) {
		
		/* 
		 * Using collections to sort ArrayList.
		 * Comparator compares fitness score of each chromosome in population.
		 * Reversed is used to sort by descending order.
		 */
		Collections.sort(population, Comparator.comparingDouble(Chromosome::getFitnessScore).reversed());
	}
	
	/**
	 * Set chromosome at specific position in population.
	 * 
	 * @param index Position in population
	 * @param chromosome Chromosome object
	 */
	public void setChromosome(int index, Chromosome chromosome) {
		this.population.add(index, chromosome);
	}
	
	/**
	 * Method to randomly shuffle chromosomes in population.
	 */
	public void shufflePopulation() {
		Collections.shuffle(this.population);
	}
	
	/**
	 * Method to calculate the fitness score of a chromosome.
	 * The lower the path distance the higher the fitness score.
	 * 
	 * @param chromosome Chromosome object
	 * @param cities Array of city objects
	 */
	public void calculateChromosomeFitness(Chromosome chromosome, double pathDistance) {
		
		/* Calculate fitness score based on path distance */
		double fitnessScore = 1 / pathDistance;
		
		/* Set chromosome fitness score */
		chromosome.setFitnessScore(fitnessScore);
	}
	
}
