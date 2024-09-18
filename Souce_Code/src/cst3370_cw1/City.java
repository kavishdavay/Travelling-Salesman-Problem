package cst3370_cw1;

/**
 * Class to store city ID and coordinates.
 * Includes a method to calculate Euclidean distance to another city.
 * 
 * @version 1.0 06 Dec 2023
 * @author M00789089
 * 
 */
public class City {
    private int cityID;
    private float xCoordinate;
    private float yCoordinate;

    
    /**
     * Initialise city object.
     * 
     * @param id City ID
     * @param x The X Coordinate
     * @param y The Y Coordinate
     */
    public City(int id, float x, float y) {
        this.cityID = id;
        this.xCoordinate = x;
        this.yCoordinate = y;
    }
    
    
    /**
     * Get the city ID.
     * 
     * @return city_ID
     */
    public int getCityID() {
        return cityID;
    }

    /**
     * Get the X coordinate.
     * 
     * @return xCoordinate
     */
    public float getXCoordinate() {
        return xCoordinate;
    }

    /**
     * Get the Y coordinate.
     * 
     * @return yCoordinate
     */
    public float getYCoordinate() {
        return yCoordinate;
    }
    
    /**
     * Calculate distance to another city.
     * 
     * @param city City object
     * @return Distance to other city
     */
    public double distanceToCity(City city) {
    	
    	/* Finding difference of X and Y coordinates between current city and city object */
    	double deltaX = this.xCoordinate - city.getXCoordinate();
    	double deltaY = this.yCoordinate - city.getYCoordinate();
    	
    	/* Applying Euclidean formula to calculate distance */
        double distance = Math.sqrt(Math.pow(deltaX, 2) + Math.pow(deltaY, 2));
        return distance;
    }
    
}
