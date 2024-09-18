package cst3370_cw1;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Genetic algorithm consisting of evaluation, selection, crossover and mutation methods.
 * It finds the best chromosome in a population based on fitness score over multiple generations.
 * Its attributes are parameters which can be fine tuned to get optimal results.
 * 
 * @version 1.0 06 Dec 2023
 * @author M00789089
 * 
 */
public class GeneticAlgorithm {

	private int maxPopulationSize;
	private double crossoverRate;
	private double mutationRate;
	private int numberOfElites;
	private int searchSpaceSize;
	
	/**
	 * Genetic algorithm constructor.
	 * 
	 * @param maxPopulationSize Population size parameter
	 * @param crossoverRate Rate at which crossover occurs
	 * @param mutationRate Rate of chromosome mutation
	 * @param numberOfElites Number of elite chromosomes in a population
	 * @param searchSpaceSize Size of search space for tournament selection
	 */
	public GeneticAlgorithm(int maxPopulationSize, double crossoverRate, double mutationRate, 
			int numberOfElites, int searchSpaceSize) {
		
		this.maxPopulationSize = maxPopulationSize;
		this.crossoverRate = crossoverRate;
		this.mutationRate = mutationRate;
		this.numberOfElites = numberOfElites;
		this.searchSpaceSize = searchSpaceSize;
	}
	
	/**
	 * Initialise the first population.
	 * 
	 * @param cities Array of city objects
	 * @return Initial Population as ArrayList of chromosome objects
	 */
	public Population startPopulation(City cities[]) {
		Population initialPopulation = new Population(this.maxPopulationSize, cities);
		return initialPopulation;
	}
	
	/**
	 * Calculates the fitness of each chromosome in a population.
	 * Sets the population's average fitness.
	 * 
	 * @param population Current population
	 * @param cities Array of city objects
	 */
	public void evaluate(Population population, City cities[]) {
    	
		for(int chromosomeIndex = 0; chromosomeIndex < population.getPopulation().size(); chromosomeIndex++) {
			Chromosome chromosome = population.getChromosome(chromosomeIndex);
			
			City[] path = searchCities(chromosome.getChromosome(), cities);	
			double pathDistance = calculatePathDistance(path);
			
			population.calculateChromosomeFitness(chromosome, pathDistance);
		}
		
		population.setAverageFitness();
		
	}

    /**
     * Checks if the number of generations parameter has been reached.
     * 
     * @param generationCount Current number of generations
     * @param maxGenerations Maximum number of generations
     * @return
     */
    public boolean isMaxGenerationsReached(int generationCount, int maxGenerations) {
        return (generationCount > maxGenerations);
    }
	
	/**
	 * Method to select second parent using tournament selection concept.
	 * Used during exploration phase.
	 * 
	 * @param population Current population
	 * @return Second parent as chromosome object
	 */
	public Chromosome tournamentSelection(Population population) {
		
		/* Create blank search space */
		Population searchSpace = new Population(this.searchSpaceSize);
		
		/* Shuffle population to randomly distribute chromosomes */
		population.shufflePopulation();
		
		/* Select chromosomes from population and add to search space */
		for(int chromosomeIndex = 0; chromosomeIndex < this.searchSpaceSize; chromosomeIndex++) {
			Chromosome candidateChromosome = population.getChromosome(chromosomeIndex);
			searchSpace.setChromosome(chromosomeIndex, candidateChromosome);
		}
		
		/* Return the chromosome with the best fitness score from search space */
		return searchSpace.getBestChromosome(0);
	}
	
	/**
	 * Method to select second parent based on roulette wheel selection.
	 * Chromosomes with a higher fitness score are more likely to be selected.
	 * Used during exploitation phase.
	 * 
	 * @param population Current population
	 * @return Chromosome object of second parent
	 */
	public Chromosome rouletteSelection(Population population) {
		
		/* Storing population chromosomes in ArrayList */
		ArrayList<Chromosome> roulettePopulation = population.getPopulation();
		
		/* 
		 * Gets population average fitness.
		 * Calculates population total fitness by multiplying with population size.
		 */
		double totalFitness = population.getAverageFitness() * roulettePopulation.size();
		
		/* 
		 * Selection threshold based on total population fitness.
		 * Random value between 0 and total fitness.
		 */
		double selectionThreshold = Math.random() * totalFitness;
		
		/* Initialise cumulative fitness */
		double cumulativeFitness = 0;
		
		/* Shuffle population randomly */
		population.shufflePopulation();
		
		/* Iterate through each chromosome in ArrayList */
		for(Chromosome chromosome : roulettePopulation) {

			/* Get chromosome fitness */
			double chromosomeFitness = chromosome.getFitnessScore();
			
			/* Add chromosome fitness to cumulative fitness */
			cumulativeFitness += chromosomeFitness;
			
			/* Return chromosome if cumulative fitness reaches selection threshold */
			if(cumulativeFitness >= selectionThreshold) {
				return chromosome;
			}
		}
		
		/* Throw error if something goes wrong */
		throw new IllegalStateException("Roulette wheel selection failed");
	}
	
	/**
	 * Crossover method to produce child chromosomes from two parents.
	 * Elitism applied whereby the best chromosomes get added directly to the next population.
	 * 
	 * @param population Current population
	 * @param selectionMethod String to determine selection method
	 * @return Population with crossover applied
	 */
	public Population crossover(Population population, String selectionMethod) {
		
		/* Initialise next population with current population size */
		Population nextPopulation = new Population(population.getPopulation().size());
		
		/* Store ArrayList of chromosomes from current population */
		ArrayList<Chromosome> chromosomes = population.getPopulation();
		
		/* Iterate through chromosomes ArrayList*/
		for(int chromosomePosition = 0; chromosomePosition < chromosomes.size(); chromosomePosition++) {
			
			/* Stores first parent chromosome */
			Chromosome firstParent = population.getBestChromosome(chromosomePosition);

			/* 
			 * Checks if crossover should be performed. 
			 * Based on crossover rate and elitism.
			 */
			if((this.crossoverRate > Math.random()) && (chromosomePosition >= this.numberOfElites)) {
				
				/* Initialise empty second parent */
				Chromosome secondParent = null;
				
				/* Generating two random points within length of first parent*/
				int point1 = (int) (Math.random() * firstParent.getChromosome().length);
				int point2 = (int) (Math.random() * firstParent.getChromosome().length);
				
				/* Getting start and end points for range of first parent genes to crossover */
				final int start = Math.min(point1, point2);
				final int end = Math.max(point1, point2);
				
				/* Checking which seletion method should be used for second parent */
				if(selectionMethod == "roulette") {
					secondParent = rouletteSelection(population);
				} else if(selectionMethod == "tournament") {
					secondParent = tournamentSelection(population);
				}
				
				/* Call orderCrossover and store chromosome returned as child */
				Chromosome child = orderCrossover(firstParent, secondParent, start, end);
				
				/* Add child chromosome to next population*/
				nextPopulation.setChromosome(chromosomePosition, child);
				
			} else {
				
				/* Add first parent directly to next population without crossover */
				nextPopulation.setChromosome(chromosomePosition, firstParent);
			}
		}
		
		/* Return next population */
		return nextPopulation;
	}
	
	/**
	 * Helper method for crossover.
	 * Returns chromosome of child produced by order crossover of two parents.
	 * 
	 * @param parent1 First parent chromosome
	 * @param parent2 Second parent chromosome
	 * @param startPoint Crossover region start index
	 * @param endPoint Crossover region end index
	 * @return Chromosome object of child
	 */
	private Chromosome orderCrossover(Chromosome parent1, Chromosome parent2, int startPoint, int endPoint) {
		
		/* Store genes of parent chromosomes */
		int[] firstParent = parent1.getChromosome();
		int[] secondParent = parent2.getChromosome();
		
		/* Initialise child genes*/
		int[] childGenes = new int[firstParent.length];
		Arrays.fill(childGenes, -1);
		
		/* Set genes from first parent into child using crossover region */
		for(int geneIndex = startPoint; geneIndex < endPoint; geneIndex++) {
			childGenes[geneIndex] = parent1.getGene(geneIndex);
		}
		
		/* Set genes from second parent into child */
		for(int geneIndex = 0; geneIndex < secondParent.length; geneIndex++) {
			
			/* Calculating gene index for second parent at end index of crossover region */
			int secondParentGeneIndex = geneIndex + endPoint;
			
			/* Circular indexing to ensure second parent gene index stays within array length */
			if(secondParentGeneIndex >= secondParent.length) {
				secondParentGeneIndex -= secondParent.length;
			}
			
			/* Select gene from second parent*/
			int chosenGene = secondParent[secondParentGeneIndex];
			
			/* Initialise match check */
			boolean matchCheck = false;

			/* Go through child genes and check if chosen gene is already in there */
			for(int gene : childGenes) {
				if(gene == chosenGene) {
					matchCheck = true;
					break;
				}
			}
			
			/* If there is no match, find empty position in child gene and add chosen gene */
			if(!matchCheck) {
				
				for(int genePosition = 0; genePosition < childGenes.length; genePosition++) {
					if(childGenes[genePosition] == -1) {
						childGenes[genePosition] = chosenGene;
						break;
					}
				}
			}
		}
		
		/* Create new child chromose with child genes */
		Chromosome child = new Chromosome(childGenes);
		
		/* Return child chromosome */
		return child;
	}
	
	/**
	 * Method to perform swap mutation on a population.
	 * 
	 * @param population Current population
	 * @return Mutated population
	 */
	public Population mutatePopulation(Population population) {
		
		/* Create blank population with size of current population */
		Population mutatedPopulation = new Population(population.getPopulation().size());
		
		/* Go through chromosomes in current population */
		for(int chromosomeIndex = 0; chromosomeIndex < population.getPopulation().size(); chromosomeIndex++) {
			
			/* Get chromosome based on best fitness scores */
			Chromosome bestChromosome = population.getBestChromosome(chromosomeIndex);
			
			/* Calculate mutation rate based on fitness score of chromosome */
			double adaptedMutationRate = adaptMutationRate(bestChromosome.getFitnessScore(), population.getAverageFitness());
			
			/* Perform only if elites have already been added to mutated population */
			if(chromosomeIndex >= this.numberOfElites) {
				
				/* Go through genes in best chromosome */
				for(int genePos = 0; genePos < bestChromosome.getChromosome().length; genePos++) {
					
					/* Check if gene mutation will occur based on adapted mutation rate */
					if(adaptedMutationRate > Math.random()) {
						
						/* Get index of gene to be swapped with */
						int swapPos = (int) (Math.random() * bestChromosome.getChromosome().length);
						
						/* Call swap mutation method */
						swapMutation(bestChromosome, genePos, swapPos);
					}
				}
			}
			
			/* Add chromosome to mutated population */
			mutatedPopulation.setChromosome(chromosomeIndex, bestChromosome);
		}
		
		/* Return mutated population */
		return mutatedPopulation;
	}
	
	/**
	 * Helper method for mutation.
	 * Dynamically adjusts mutation rate based on the fitness score of chromosome.
	 * 
	 * @param fitnessScore Chromosome fitness score
	 * @param averageFitness Average fitness of current population
	 * @return
	 */
	private double adaptMutationRate(double fitnessScore, double averageFitness) {
		
		/* Check if fitness score is greater than average and is not zero */
		if((fitnessScore > averageFitness) && (averageFitness != 0)) {
			
			/* Initialise adapted rate */
			double adaptedRate;
			
			/* Find difference ratio between fitness score and average fitness */
			double fitnessDifferenceRatio = (fitnessScore - averageFitness) / averageFitness;
			
			/* Sigmoid function to get value between 0 and 1 */
			double scaledRatio = 1 / (1 + Math.exp(-fitnessDifferenceRatio));
			
			/* Calculate adapted rate by multiplying mutation rate with ratio */
			adaptedRate = this.mutationRate * scaledRatio;
			return adaptedRate;
		} else {
			
			/* Mutation rate for chromosomes with lower than average fitness is not adjusted */
			return this.mutationRate;
		}
	}
	
	/**
	 * Helper method to swap genes in a chromosome.
	 * 
	 * @param chromosome Chromosome object
	 * @param genePosition First gene position
	 * @param swapPosition Second gene position
	 */
	private void swapMutation(Chromosome chromosome, int genePosition, int swapPosition) {
		
		/* Initialise genes to be swapped */
		int firstGene = chromosome.getGene(genePosition);
		int secondGene = chromosome.getGene(swapPosition);
		
		/* Swap genes */
		chromosome.setGene(swapPosition, firstGene);
		chromosome.setGene(genePosition, secondGene);
	}
	
	/**
	 * 2 Opt local search to improve chromosomes with the best fitness scores.
	 * 
	 * @param population Current population
	 * @param cities Array of city objects
	 * @param twoOptThreshold Ratio of population to perform 2opt on
	 * @return Population of chromosomes with 2opt applied
	 */
	public Population performTwoOptLocalSearch(Population population, City cities[], double twoOptThreshold) {
		
		/* Creating new population to store optimised chromosomes */
		Population optimisedPopulation = new Population(population.getPopulation().size());
		
		/* Calculating 2 opt point from threshold and current population size */
		int populationSize = population.getPopulation().size();
		int twoOptPoint = (int) Math.round(populationSize * twoOptThreshold);
		
		/* Getting chromosomes in current population to perform 2opt */
		for(int chromosomeIndex = 0; chromosomeIndex < twoOptPoint; chromosomeIndex++) {
			
			/* Sort population by best chromosome fitness scores and get chromosome at index */
			Chromosome currentChromosome = population.getBestChromosome(chromosomeIndex);	
			
			optimiseChromosomeTwoOpt(population, currentChromosome, cities);
			optimisedPopulation.setChromosome(chromosomeIndex, currentChromosome);
		}
		
		/* Number of chromosomes that will not have 2opt performed */
		int restOfPopulation = populationSize - twoOptPoint;
		
		/* Add the rest of chromosomes from current population without performing 2opt */
	    for (int chromosomeIndex = 0; chromosomeIndex < restOfPopulation; chromosomeIndex++) {
	        optimisedPopulation.setChromosome(chromosomeIndex + twoOptPoint, population.getChromosome(chromosomeIndex));
	    }
		
		return optimisedPopulation;
	}
	
	/**
	 * Helper method to select range of genes for 2opt.
	 * Changes gene in current chromosome if path distance decreases after reversal.
	 * 
	 * @param population
	 * @param currentChromosome
	 * @param cities
	 */
	private void optimiseChromosomeTwoOpt(Population population, Chromosome currentChromosome, City cities[]) {
		
		/* Initialise variables with chromosome information */
		int[] currentGenes = currentChromosome.getChromosome();
		double currentDistance = currentChromosome.getFitnessScore() / 1;
		
		/* Going through edges in chromosome and performing 2 opt */
		for(int edge1 = 1; edge1 < (currentGenes.length); edge1++) {
			for(int edge2 = (edge1 + 1); edge2 < currentGenes.length; edge2++) {
				
				/* Gets path after 2opt */
				int[] genesTwoOpt = applyTwoOpt(currentGenes, edge1, edge2);
				
				/* Calculating distance for new path */
				City[] pathTwoOpt = searchCities(genesTwoOpt, cities);
				double distanceTwoOpt = calculatePathDistance(pathTwoOpt);
				
				/* Genes after 2opt set to current chromosome if new path distance is improved */
				if(distanceTwoOpt < currentDistance) {
					for(int geneIndex = 0; geneIndex < genesTwoOpt.length; geneIndex++) {
						currentChromosome.setGene(geneIndex, genesTwoOpt[geneIndex]);
					}
				}

			}
		}
	}
	
	/**
	 * Helper method to perform reversal on selected range for 2opt.
	 * 
	 * @param genes Array of chromosome genes
	 * @param edge1 Start index of reversal range
	 * @param edge2 End index of reversal range
	 * @return
	 */
	private int[] applyTwoOpt(int[] genes, int edge1, int edge2) {
		
		/* Initialise empty array */
		int[] genesTwoOpt = new int[genes.length];
		
		/* Adding genes before start edge to two opt array */
		for(int geneIndex = 0; geneIndex < edge1; geneIndex++) {
			genesTwoOpt[geneIndex] = genes[geneIndex];
		}
		
	    /* Reversing genes between edges and adding it to array */
	    for (int geneIndex = edge2, newIndex = edge1; geneIndex >= edge1; geneIndex--, newIndex++) {
	        genesTwoOpt[newIndex] = genes[geneIndex];
	    }
		
		/* Adding genes after end edge to two opt array */
		for(int geneIndex = (edge2 + 1); geneIndex < genes.length; geneIndex++) {
			genesTwoOpt[geneIndex] = genes[geneIndex];
		}

		return genesTwoOpt;
	}
	
	/**
	 * Finds corresponding city objects based on gene from chromosome.
	 * 
	 * @param chromosome Chromosome object
	 * @param cities Array of city objects
	 * @return Corresponding array of city objects
	 */
	private City[] searchCities(int[] genes, City cities[]) {
		
		/* Storing chromosome genes in array */
		int[] pathArray = genes;
		
		/* Initialise new array of city objects with length of chromosome */
		City[] citiesArray = new City[pathArray.length];
		
		/* Counter for number of matching city IDs*/
		int matchCount = 0;
		
		/* Go through each gene in the chromosome to find match in cities array */
		for(int geneIndex = 0; geneIndex < pathArray.length; geneIndex++) {
			int currentGene = pathArray[geneIndex];		// Get current gene
			boolean matchFound = false;						// Initialise check
			
			/* Go through cities array until a match is found */
			for(int cityIndex = 0; (cityIndex < cities.length) && (matchFound == false); cityIndex++) {
				
				/* Get the city ID */
				int cityID = cities[cityIndex].getCityID();
				
				/* If the cityID corresponds to the gene, add the city to the citiesArray */
				if(currentGene == cityID) {
					citiesArray[matchCount] = cities[cityIndex];
					matchCount++;
					matchFound = true;						// Set check to true when found
				}
			}
		}
		
		return citiesArray;
	}

	/**
	 * Method to calculate the total distance of a path.
	 * 
	 * @param path Array of city objects
	 * @return double totalDistance
	 */
	private double calculatePathDistance(City path[]) {
		
		/* Initialise total distance */
		double totalDistance = 0;
		
		/* Go through each city in path, calculate distance to next city and add to total distance */
		for(int cityIndex = 0; cityIndex < path.length - 1; cityIndex++) {
			
			City currentCity = path[cityIndex];
			City nextCity = path[cityIndex + 1];
			
			double distance = currentCity.distanceToCity(nextCity);
			
			totalDistance += distance;
		}
		
		/* Complete cycle by adding distance from last city to first city to total distance */
		totalDistance += path[path.length - 1].distanceToCity(path[0]);		
		
		return totalDistance;
	}
	
}