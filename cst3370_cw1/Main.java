package cst3370_cw1;

import java.util.Scanner;

/**
 * The main source code for the Travelling Salesman problem.
 * Genetic algorithm with 2opt local search is used to find the optimal travel path.
 * Gets cities from a file.
 * 
 * @version 1.0 06 Dec 2023
 * @author M00789089
 * 
 */
public class Main {
	
	/* Parameters */
	static final int MAX_POPULATION_SIZE = 150;
	static final double CROSSOVER_RATE = 0.9;
	static final double MUTATION_RATE = 0.001;
	static final int NUMBER_OF_ELITES = 4;
	static final int SEARCH_SPACE_SIZE = 10;
	static final int MAX_GENERATIONS = 800;
	static final double EXPLORATION_THRESHOLD = 0.9;
	static final double TWO_OPT_THRESHOLD = 0.1;
	
	
	public static void main(String[] args) {
		
		/* Scanner to get user input */
		Scanner scanner = new Scanner(System.in);
		
		/* Initialise variables */
		String sampleFile = null;
		int[][] readOutput = null;
		
		/* Asks for file until file is found and in correct format */
        boolean foundCheck = false;
        while (!foundCheck) {
            System.out.print("Enter file path: ");
            sampleFile = scanner.nextLine();

            try {
            	
            	/* Handler to read files */
                readOutput = FileHandler.read(sampleFile);
                foundCheck = true;
            } catch (RuntimeException e) {
            	
            	/* Prints error if exception occurs */
                System.out.println(e.getMessage() + "\n\n");
            }
        }
        
        /* Closing scanner */
        scanner.close();
        
        /* Initializing cities from file output */
		City[] cities = initializeCities(readOutput);
		
		/* Start timer */
		final long startTime = System.nanoTime();
		
		/* Initialise genetic algorithm with parameters */
		GeneticAlgorithm geneticAlgorithm = new GeneticAlgorithm(MAX_POPULATION_SIZE, CROSSOVER_RATE, MUTATION_RATE, 
				NUMBER_OF_ELITES, SEARCH_SPACE_SIZE);
		
		/* Initialise start population */
		Population population = geneticAlgorithm.startPopulation(cities);
		
		/* Evaluate start population */
		geneticAlgorithm.evaluate(population, cities);
		
		/* Calculate best distance in start population */
		double startBestFitness = population.getBestChromosome(0).getFitnessScore();
		double startDistance = 1 / startBestFitness;
		
		/* Store start population best path */
		int [] startPath = population.getBestChromosome(0).getChromosome();
		
		/* Calling method to print start population best path and distance */
		printInitialPath(startPath, startDistance);
		
		/* Counter for number of generations */
		int generationCount = 1;
		
		/* Loop to run genetic algorithm until maximum generations */
		while(!geneticAlgorithm.isMaxGenerationsReached(generationCount, MAX_GENERATIONS)) {
			
			/* Apply 2 opt local search to population */
			population = geneticAlgorithm.performTwoOptLocalSearch(population, cities, TWO_OPT_THRESHOLD);
			
			/* Check which selection method should be used for crossover based on number of generations */
			if(generationCount > (MAX_GENERATIONS * EXPLORATION_THRESHOLD)) {
				
				/* Exploitation stage */
				population = geneticAlgorithm.crossover(population, "roulette");
			} else {
				
				/* Exploration stage*/
				population = geneticAlgorithm.crossover(population, "tournament");
			}
			
			/* Apply mutation to population */
			population = geneticAlgorithm.mutatePopulation(population);
			
			geneticAlgorithm.evaluate(population, cities);
			
			/* Increment generation counter */
			generationCount++;
		}
		
		/* Get current time and calculate time taken in nanoseconds */
		long currentTime = System.nanoTime();
		double timeTaken = (currentTime - startTime);
		
		/* Get best path in final population */
		Chromosome bestChromosome = population.getBestChromosome(0);
		int [] bestPath = bestChromosome.getChromosome();
		
		/* Calculate path distance */
		double bestFitness = bestChromosome.getFitnessScore();
		double bestDistance = 1 / bestFitness;
		
		/* Print best path */
		printOptimalPath(bestPath, bestDistance);
		
		/* Print time taken */
		System.out.println("Time taken: " + timeTaken + " Nanoseconds");
	}
	
    /**
     * Method to initialise city objects from 2D array.
     * 
     * @param readOutput 2D array from file
     * @return Array of city objects
     */
    private static City[] initializeCities(int [][] readOutput) {
        
    	int numberOfCities = readOutput.length;
        City[] cities = new City[numberOfCities];
        
        for (int i = 0; i < numberOfCities; i++) {
            cities[i] = new City(readOutput[i][0], readOutput[i][1], readOutput[i][2]);
        }
        
        return cities;
    }
    
    /**
     * Helper method to print the initial path and distance.
     * 
     * @param startPath Array of city IDs
     * @param startDistance Distance of path
     */
    private static void printInitialPath(int [] startPath, double startDistance) {

    	System.out.println("Best initial Path: ");
        for (int cityIndex = 0; cityIndex < startPath.length; cityIndex++) {
            System.out.print(startPath[cityIndex] + " -> ");
        }
        System.out.print(startPath[0]);
		System.out.println("\n\nInitial path distance = " + startDistance + "\n\n");

    }
    
    /**
     * Helper method to print optimal path and distance.
     * 
     * @param bestPath Array of city IDs
     * @param bestDistance Distance of path
     */
    private static void printOptimalPath(int [] bestPath, double bestDistance) {
        System.out.println("Best Path: ");
        for (int cityIndex = 0; cityIndex < bestPath.length; cityIndex++) {
            System.out.print(bestPath[cityIndex] + " -> ");
        }
        System.out.print(bestPath[0]);
        
        System.out.println();
        System.out.println("\nBest Distance: " + bestDistance);
    }
}