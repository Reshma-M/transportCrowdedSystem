/*
 * System configuration:
 * 1. Store the compartment information
 *
 *
 */


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <map>
#include <cmath>
#include <stdio.h>
#include <time.h>

using namespace std;

const long NUMBER_OF_UPDATES_TO_SETTLE = 100000000;
const long NUMBER_OF_UPDATES = NUMBER_OF_UPDATES_TO_SETTLE * 1;
int NUMBER_OF_SNAPSHOTS = 1500;
const int NUMBER_OF_VESICLES = 200;
const int NUMBER_OF_SITES = 1000;
const int NUMBER_OF_MICROTUBULES = 10;
const int NUMBER_OF_TYPES = 3;
const int NUMBER_OF_EXCLUDED_NEIGHBOURS = 3;
const int LENGTH_OF_BREAK = 13;
const double SNAPSHOT_INTERVAL = 0.2;
const double STND_PROPENSITY = 0;
const double Sm_HOPPING_PROPENSITY = 100;
const double StA_HOPPING_PROPENSITY = 40;
const double SmA_CHANGETRACK_PROPENSITY = 0.2/10;
const double StA_CHANGETRACK_PROPENSITY = 0.2/10;
const double SmR_CHANGETRACK_PROPENSITY = 0.2/10;
const double CHANGETYPE_FACTOR = 1;
const double Sm_To_StA_CHANGETYPE_PROPENSITY = 0.08 * CHANGETYPE_FACTOR;
const double StA_To_SmA_CHANGETYPE_PROPENSITY = 0.8;
const double StA_To_SmR_CHANGETYPE_PROPENSITY = 0.64;
const int TURNBACK_FACTOR = 1;
const int SIDESTEP_FACTOR = 1;

//const double Sm_To_Sm_CHANGETYPE_PROPENSITY = 0.1; // To test the cases where there are only Smooth Vesicles;

/*Initializations for Metadata*/

bool IS_SYSTEM_SETTLED = false;
long THIS_UPDATE = 0;
double TIME_SINCE_SNAPSHOT = 0;
double TOTAL_TIME = 0;
long NUMBER_OF_HOPS = 0;
long NUMBER_OF_CHANGETRACKS = 0;
long NUMBER_OF_CHANGETYPES = 0;
vector<long> CURRENT_AT_LOCATION_1;
vector<long> CURRENT_AT_LOCATION_2;
vector<long> CURRENT_AT_LOCATION_3;
vector<vector<long> > DENSITY;
vector<vector<long> > MICROTUBULE_CONTINUING_CLUSTERS;
vector<long> CONTINUING_CLUSTERS;
map<long, long> CLUSTER_DURATIONS;
map<long, long> CHANGETRACK_FREQUENCY;
map<long, long> CHANGETYPE_FREQUENCY;
vector<map<long, long> > MICROTUBULE_CLUSTER_DURATIONS;

long SmA2StAAtSC = 0;
long StA2SmAAtSC = 0;
long StA2SmRAtSC = 0;
long SmR2StAAtSC = 0;
long SmA2StAAwayFromSC = 0;
long StA2SmAAwayFromSC = 0;
long StA2SmRAwayFromSC = 0;
long SmR2StAAwayFromSC = 0;
long AllStA2SmA = 0;
long AllStA2SmR = 0;
long AllSmA2StA = 0;
long AllSmR2StA = 0;

/*
 * Enum of the directions of Vesicle motion
 */
enum Direction
{
    ANTEROGRADE = -1, RETROGRADE = 1
};

/*
 * Enum of the neighbours of a microtubule
 */
enum NeighbourMicrotubule
{
    PREVIOUS = -1, NEXT = 1
};

/*
 * Enum of Vesicle Types
 */
enum VesicleType
{
    SMOOTH_ANTEROGRADE, STAGGERED_ANTEROGRADE, SMOOTH_RETROGRADE, STANDSTILL
};

/*
 * Class that models a Vesicle. vesicleId uniquely identifies Vesicles. Values
 * of vesicleId start at 1.
 */
class Vesicle
{
    long vesicleId;
    int microtubuleId;
    int siteId;
    Direction direction;
    VesicleType vesicleType;
    long numberOfChangetypes;
    long numberOfChangetracks;
public:

    long getVesicleId()
    {
        return (vesicleId);
    }

    void setVesicleId(long vId)
    {
        vesicleId = vId;
    }

    int getMicrotubuleId()
    {
        return (microtubuleId);
    }

    void setMicrotubuleId(int mId)
    {
        microtubuleId = mId;
    }

    int getSiteId()
    {
        return (siteId);
    }

    void setSiteId(int sId)
    {
        siteId = sId;
    }

    Direction getDirection()
    {
        return (direction);
    }

    void setDirection()
    {
        direction = ANTEROGRADE;
        if (getType() == SMOOTH_RETROGRADE)
        {
            direction = RETROGRADE;
        }
    }

    VesicleType getType()
    {
        return (vesicleType);
    }

    void setType(VesicleType vT)
    {
        vesicleType = vT;
        setDirection();
    }

    long getNumberOfChangetypes()
    {
        return (numberOfChangetypes);
    }

    void setNumberOfChangetypes(long cty)
    {
        numberOfChangetypes = cty;
    }

    long getNumberOfChangetracks()
    {
        return (numberOfChangetracks);
    }

    void setNumberOfChangetracks(long ctr)
    {
        numberOfChangetracks = ctr;
    }

    Vesicle(long vId, int mId, int sId, VesicleType vesicleType, long cty, long ctr)
    {
        setVesicleId(vId);
        setMicrotubuleId(mId);
        setSiteId(sId);
        setType(vesicleType);
        setNumberOfChangetypes(cty);
        setNumberOfChangetracks(ctr);
    }
};

/*
 * Needed for the set container
 */
class Comp
{
public:

    bool operator()(Vesicle v1, Vesicle v2)
    {
        if (v1.getMicrotubuleId() < v2.getMicrotubuleId())
        {
            return (true);
        }
        else
        {
            if (v1.getMicrotubuleId() == v2.getMicrotubuleId())
            {
                if (v1.getSiteId() < v2.getSiteId())
                {
                    return (true);
                }
            }
            return (false);
        }
    }
};

set<Vesicle, Comp> VESICLES;

/*
 * Data Structure to hold the location of a Vesicle. Helps to return both the Ids -
 * microtubule and site - from a function
 */
struct Location
{
    int microtubuleId;
    int siteId;
};

/*
 * Returns a random number between 0 and 1 (including 0)
 */
double randomNumber()
{
    return (rand() / (double) RAND_MAX);
}

/*
 * Returns the next Vesicle Id. The first Vesicle gets the id 1
 */
long generateVesicleId()
{
    static long initialId = 0;
    return (++initialId);
}

/*
 * True if the Site is occupied; false otherwise
 */
bool isSiteOccupied(Location location)
{
    Vesicle vesicle = Vesicle(0, 0, 0, SMOOTH_ANTEROGRADE, 0, 0);
    vesicle.setMicrotubuleId(location.microtubuleId);
    vesicle.setSiteId(location.siteId);
    if (VESICLES.find(vesicle) != VESICLES.end())
    {
        return (true);
    }
    else
    {
        return (false);
    }
}

/*
 * False if any of the sites in listOfSites is occupied; true otherwise.
 */
bool isListofSitesEmpty(vector<Location> listOfSites)
{
    for (unsigned int i = 0; i < listOfSites.size(); i++)
    {
        if (isSiteOccupied(listOfSites.at(i)))
        {
            return (false);
        }
    }
    return (true);
}

/*
 * True if a Vesicle can hop into the given site; false otherwise
 */
bool isSiteHoppable(int microtubuleNumber, int siteNumber, Direction direction)
{
    int neighboursToBeChecked = NUMBER_OF_EXCLUDED_NEIGHBOURS;
    vector<Location> sitesToBeChecked;
    while (neighboursToBeChecked >= 0)
    {
        Location location;
        location.microtubuleId = microtubuleNumber;
        location.siteId = (siteNumber + NUMBER_OF_SITES + (direction * neighboursToBeChecked)) % NUMBER_OF_SITES;
        sitesToBeChecked.push_back(location);
        neighboursToBeChecked--;
    }
    return (isListofSitesEmpty(sitesToBeChecked));
}

/*
 * Checks if a Vesicle can be placed at the given site.
 */
bool isSiteOccupiable(int microtubuleNumber, int siteNumber)
{
    return (isSiteHoppable(microtubuleNumber, siteNumber, ANTEROGRADE) && isSiteHoppable(microtubuleNumber, siteNumber, RETROGRADE));
}

/*
 * Updates the siteNumber of the vesicle and places it in the container.
 */
void placeVesicle(Vesicle v, int siteNumber)
{
    v.setSiteId(siteNumber);
    VESICLES.insert(v);
}

/*
 * Fill the Axon with Vesicles so that they will not be able to hop or changeTrack
 */
/*void stuffVesicles()
 {
 for (int i = 0; i < NUMBER_OF_MICROTUBULES; i++)
 {
 for (int j = 0; j < NUMBER_OF_SITES; j++)
 {
 if (isSiteOccupiable(i, j) == true)
 {
 VESICLES.insert(Vesicle(generateVesicleId(), i, j, SMOOTH_ANTEROGRADE));
 }
 }
 }
 cout << "The number of Vesicles is " << VESICLES.size() << endl;
 }*/

/*
 * Creates NUMBER_OF_VESICLES number of Vesicles. They are placed at random locations.
 */
void makeVesicles()
{
    int numberOfVesicles = NUMBER_OF_VESICLES;
    while (numberOfVesicles > 0)
    {
        int microtubuleNumber = rand() % NUMBER_OF_MICROTUBULES;
        int siteNumber = rand() % NUMBER_OF_SITES;
        VesicleType vT = SMOOTH_ANTEROGRADE;
        if (numberOfVesicles <= (NUMBER_OF_VESICLES / 2))
        {
            vT = SMOOTH_RETROGRADE;
        }
        if (numberOfVesicles <= NUMBER_OF_VESICLES / 10)
        {
            vT = STAGGERED_ANTEROGRADE;
        }
        if (isSiteOccupiable(microtubuleNumber, siteNumber))
        {
            VESICLES.insert(Vesicle(generateVesicleId(), microtubuleNumber, siteNumber, vT, 0, 0));
            numberOfVesicles--;
        }
    }
}

/*
 * Create breaks in the given locations
 */
void makeBreak(int microtubuleNumber, int startSite)
{
    for (int i = 0; i < LENGTH_OF_BREAK; i++)
    {
        VESICLES.insert(Vesicle(generateVesicleId(), microtubuleNumber, (startSite + i) % NUMBER_OF_SITES, STANDSTILL, 0, 0));
    }
}

/*
 * Returns the propensity for the Vesicle to change type
 */
double getChangeTypePropensity(Vesicle v)
{
    if (v.getType() == STANDSTILL)
    {
        return (0);
    }
    if (v.getType() == STAGGERED_ANTEROGRADE)
    {
        return (StA_To_SmA_CHANGETYPE_PROPENSITY + StA_To_SmR_CHANGETYPE_PROPENSITY);
    }
    else
    {
        return (Sm_To_StA_CHANGETYPE_PROPENSITY);
        //return (Sm_To_Sm_CHANGETYPE_PROPENSITY);
    }
}

/*
 * Returns the Side-Stepping propensity of the given Vesicle
 */
double getChangeTrackPropensity(Vesicle v, NeighbourMicrotubule neighbour)
{
    if (v.getType() == STANDSTILL)
    {
        return (0);
    }
    if (isSiteOccupiable((v.getMicrotubuleId() + NUMBER_OF_MICROTUBULES + (1 * neighbour)) % NUMBER_OF_MICROTUBULES, v.getSiteId()))
    {
        if (v.getType() == SMOOTH_ANTEROGRADE)
        {
            return (SmA_CHANGETRACK_PROPENSITY);
        }
        if (v.getType() == STAGGERED_ANTEROGRADE)
        {
            return (StA_CHANGETRACK_PROPENSITY);
        }
        if (v.getType() == SMOOTH_RETROGRADE)
        {
            return (SmR_CHANGETRACK_PROPENSITY);
        }
    }
    return (0);
}

/*
 * Returns the hopping propensity of the given Vesicle
 */
double getHoppingPropensity(Vesicle v)
{
    if (v.getType() == STANDSTILL)
    {
        return (0);
    }
    if (isSiteHoppable(v.getMicrotubuleId(), (v.getSiteId() + NUMBER_OF_SITES + (1 * v.getDirection())) % NUMBER_OF_SITES, v.getDirection()))
    {
        if (v.getType() == STAGGERED_ANTEROGRADE)
        {
            return (StA_HOPPING_PROPENSITY);
        }
        else
        {
            return (Sm_HOPPING_PROPENSITY);
        }
    }
    return (0);
}

/*
 * Returns the propensity for the given Vesicle
 */
double getVesiclePropensity(Vesicle v)
{
    if (getHoppingPropensity(v) > 0)
    {
        return (getHoppingPropensity(v) + (getChangeTypePropensity(v) / TURNBACK_FACTOR) + (getChangeTrackPropensity(v, PREVIOUS) / SIDESTEP_FACTOR) + (getChangeTrackPropensity(v, NEXT) / SIDESTEP_FACTOR));
    }
    else
    {
        double propensity = 0;
        propensity += getChangeTrackPropensity(v, PREVIOUS);
        propensity += getChangeTrackPropensity(v, NEXT);
        propensity += getChangeTypePropensity(v);
        return (propensity);
    }
}

/*
 * Updates current
 */
void updateCurrent(Vesicle v)
{
    int siteNumber = v.getSiteId();
    if (v.getType() == SMOOTH_RETROGRADE)
    {
        if (siteNumber == (NUMBER_OF_SITES / 4))
        {
            CURRENT_AT_LOCATION_1.at(v.getType()) += 1;
        }
        if (siteNumber == (NUMBER_OF_SITES / 2))
        {
            CURRENT_AT_LOCATION_2.at(v.getType()) += 1;
        }
        if (siteNumber == ((NUMBER_OF_SITES * 3) / 4))
        {
            CURRENT_AT_LOCATION_3.at(v.getType()) += 1;
        }
    }
    else
    {
        if (siteNumber == ((NUMBER_OF_SITES / 4) + 1))
        {
            CURRENT_AT_LOCATION_1.at(v.getType()) += 1;
        }
        if (siteNumber == ((NUMBER_OF_SITES / 2) + 1))
        {
            CURRENT_AT_LOCATION_2.at(v.getType()) += 1;
        }
        if (siteNumber == (((NUMBER_OF_SITES * 3) / 4) + 1))
        {
            CURRENT_AT_LOCATION_3.at(v.getType()) += 1;
        }
    }
}

/*
 * Updates the density
 */
void updateDensity()
{
    if (THIS_UPDATE % NUMBER_OF_VESICLES == 0)
    {
        for (set<Vesicle>::iterator it = VESICLES.begin(); it != VESICLES.end(); it++)
        {
            Vesicle v = (*it);
            if (v.getType() != STANDSTILL)
            {
                DENSITY.at(v.getType()).at(v.getSiteId()) += 1;
            }
        }
    }
}

/*
 * This function takes a snapshot of the axon. Each microtubule is written to a separate file.
 * In the snapshot, the content of a site is the vesicleId of the vesicle in that site, if
 * that site is occupied or -1 otherwise.
 */
void writeSnapshot()
{
    for (int i = 0; i < NUMBER_OF_MICROTUBULES; i++)
    {
        ostringstream snapshot;
        ostringstream typeSnapshot;
        for (int j = 0; j < NUMBER_OF_SITES; j++)
        {
            Location location;
            location.microtubuleId = i;
            location.siteId = j;
            if (isSiteOccupied(location))
            {
                Vesicle vesicle = Vesicle(0, 0, 0, SMOOTH_ANTEROGRADE, 0, 0);
                vesicle.setMicrotubuleId(location.microtubuleId);
                vesicle.setSiteId(location.siteId);
                set<Vesicle>::iterator it = VESICLES.find(vesicle);
                Vesicle v = (*it);
                if (v.getType() != STANDSTILL)
                {
                    snapshot << v.getVesicleId() << ",";
                    typeSnapshot << v.getType() << ",";
                }
                else
                {
                    snapshot << "-1,";
                    typeSnapshot << "-1,";
                }
            }
            else
            {
                snapshot << "-1,";
                typeSnapshot << "-1,";
            }
        }
        ostringstream stringStream;
        ostringstream typeStringStream;
        stringStream << "kymo" << i << ".csv";
        typeStringStream << "kymoType" << i << ".csv";
        string fileName;
        string typeFileName;
        fileName = stringStream.str();
        typeFileName = typeStringStream.str();
        ofstream microtubuleWriter;
        ofstream microtubuleTypeWriter;
        microtubuleWriter.open(fileName.c_str(), ios::out | ios::app);
        microtubuleTypeWriter.open(typeFileName.c_str(), ios::out | ios::app);
        microtubuleWriter << snapshot.str() << endl;
        microtubuleTypeWriter << typeSnapshot.str() << endl;
        microtubuleWriter.close();
        microtubuleTypeWriter.close();
    }
}

/*
 * Takes a 2d snapshot of the axon
 */
vector<long> takeClusterSnapshot()
{
    vector<long> currentClusters;
    currentClusters.assign(NUMBER_OF_SITES, 0);
    for (set<Vesicle>::iterator it = VESICLES.begin(); it != VESICLES.end(); it++)
    {
        Vesicle v = (*it);
        if (v.getType() != STANDSTILL)
        {
            currentClusters.at(v.getSiteId()) += 1;
        }
    }
    return (currentClusters);
}

/*
 * Takes the snapshot of the microtubule passed as parameter
 */
vector<long> takeMicrotubuleClusterSnapshot(int mId)
{
    vector<long> currentClusters;
    currentClusters.assign(NUMBER_OF_SITES, 0);
    for (set<Vesicle>::iterator it = VESICLES.begin(); it != VESICLES.end(); it++)
    {
        Vesicle v = (*it);
        if (v.getType() != STANDSTILL && v.getMicrotubuleId() == mId)
        {
            currentClusters.at(v.getSiteId()) += 1;
        }
    }
    return (currentClusters);
}

/*
 * Updates the clusters in each microtubule
 */

void updateMicrotubuleClusters()
{
    for (int j = 0; j < NUMBER_OF_MICROTUBULES; j++)
    {
        vector<long> clusterSnapshot = takeMicrotubuleClusterSnapshot(j);
        for (int i = 0; i < NUMBER_OF_SITES; i++)
        {
            if ((MICROTUBULE_CONTINUING_CLUSTERS.at(j).at(i) >= -1) && (clusterSnapshot.at(i) > 0))
            {
                MICROTUBULE_CONTINUING_CLUSTERS.at(j).at(i) += 1;
            }
            if ((MICROTUBULE_CONTINUING_CLUSTERS.at(j).at(i) > -1) && (clusterSnapshot.at(i) == 0))
            {
                if (MICROTUBULE_CLUSTER_DURATIONS.at(j).count(MICROTUBULE_CONTINUING_CLUSTERS.at(j).at(i)) == 0)
                {
                    MICROTUBULE_CLUSTER_DURATIONS.at(j)[MICROTUBULE_CONTINUING_CLUSTERS.at(j).at(i)] = 1;
                }
                else
                {
                    MICROTUBULE_CLUSTER_DURATIONS.at(j).at(MICROTUBULE_CONTINUING_CLUSTERS.at(j).at(i)) += 1;
                }
                MICROTUBULE_CONTINUING_CLUSTERS.at(j).at(i) = -1;
            }
        }
    }
}

/*
 * Updates Clusters
 */
void updateClusters()
{
    if (TIME_SINCE_SNAPSHOT > SNAPSHOT_INTERVAL)
    {
        vector<long> clusterSnapshot = takeClusterSnapshot();
        for (int i = 0; i < NUMBER_OF_SITES; i++)
        {
            if ((CONTINUING_CLUSTERS.at(i) >= -1) && (clusterSnapshot.at(i) > 0))
            {
                CONTINUING_CLUSTERS.at(i) += 1;
            }
            if ((CONTINUING_CLUSTERS.at(i) > -1) && (clusterSnapshot.at(i) == 0))
            {
                if (CLUSTER_DURATIONS.count(CONTINUING_CLUSTERS.at(i)) == 0)
                {
                    CLUSTER_DURATIONS[CONTINUING_CLUSTERS.at(i)] = 1;
                }
                else
                {
                    CLUSTER_DURATIONS.at(CONTINUING_CLUSTERS.at(i)) += 1;
                }
                CONTINUING_CLUSTERS.at(i) = -1;
            }
        }
        TIME_SINCE_SNAPSHOT = 0;
        updateMicrotubuleClusters();
    }
}

void updateProductiveChangeTypeCount(Vesicle oldVesicle, Vesicle newVesicle)
{
    if (getHoppingPropensity(newVesicle) == 0)
    {
        return;
    }
    int siteId = oldVesicle.getSiteId();
    VesicleType oldType = oldVesicle.getType();
    VesicleType newType = newVesicle.getType();
    if (siteId > 400 && siteId < 600)
    {
        if (oldType == SMOOTH_ANTEROGRADE && newType == STAGGERED_ANTEROGRADE)
        {
            SmA2StAAtSC++;
        }
        if (oldType == SMOOTH_RETROGRADE && newType == STAGGERED_ANTEROGRADE)
        {
            SmR2StAAtSC++;
        }
        if (oldType == STAGGERED_ANTEROGRADE && newType == SMOOTH_ANTEROGRADE)
        {
            StA2SmAAtSC++;
        }
        if (oldType == STAGGERED_ANTEROGRADE && newType == SMOOTH_RETROGRADE)
        {
            StA2SmRAtSC++;
        }
    }
    if (siteId > 50 && siteId < 250)
    {
        if (oldType == SMOOTH_ANTEROGRADE && newType == STAGGERED_ANTEROGRADE)
        {
            SmA2StAAwayFromSC++;
        }
        if (oldType == SMOOTH_RETROGRADE && newType == STAGGERED_ANTEROGRADE)
        {
            SmR2StAAwayFromSC++;
        }
        if (oldType == STAGGERED_ANTEROGRADE && newType == SMOOTH_ANTEROGRADE)
        {
            StA2SmAAwayFromSC++;
        }
        if (oldType == STAGGERED_ANTEROGRADE && newType == SMOOTH_RETROGRADE)
        {
            StA2SmRAwayFromSC++;
        }
    }
    if (oldType == SMOOTH_ANTEROGRADE && newType == STAGGERED_ANTEROGRADE)
    {
        AllSmA2StA++;
    }
    if (oldType == SMOOTH_RETROGRADE && newType == STAGGERED_ANTEROGRADE)
    {
        AllSmR2StA++;
    }
    if (oldType == STAGGERED_ANTEROGRADE && newType == SMOOTH_ANTEROGRADE)
    {
        AllStA2SmA++;
    }
    if (oldType == STAGGERED_ANTEROGRADE && newType == SMOOTH_RETROGRADE)
    {
        AllStA2SmR++;
    }
}

/*
 * Remove the currentVesicle from the set and add updatedVesicle
 */
void updateVesicle(Vesicle currentVesicle, Vesicle updatedVesicle)
{
    if (currentVesicle.getType() == STANDSTILL)
    {
        cout << "This vesicle should not be updated" << endl;
    }
    if (IS_SYSTEM_SETTLED == true)
    {
        updateDensity();
        updateClusters();
        updateProductiveChangeTypeCount(currentVesicle, updatedVesicle);
        // The following check is necessary so that changeTypes and changeTracks dont result in the current being updated
        if ((currentVesicle.getMicrotubuleId() == updatedVesicle.getMicrotubuleId()) && (currentVesicle.getType() == updatedVesicle.getType()))
        {
            updateCurrent(updatedVesicle);
            if (currentVesicle.getVesicleId() != updatedVesicle.getVesicleId())
            {
                if (CHANGETYPE_FREQUENCY.count(currentVesicle.getNumberOfChangetypes()) == 0)
                {
                    CHANGETYPE_FREQUENCY[currentVesicle.getNumberOfChangetypes()] = 1;
                }
                else
                {
                    CHANGETYPE_FREQUENCY.at(currentVesicle.getNumberOfChangetypes()) += 1;
                }
                if (CHANGETRACK_FREQUENCY.count(currentVesicle.getNumberOfChangetracks()) == 0)
                {
                    CHANGETRACK_FREQUENCY[currentVesicle.getNumberOfChangetracks()] = 1;
                }
                else
                {
                    CHANGETRACK_FREQUENCY.at(currentVesicle.getNumberOfChangetracks()) += 1;
                }
                updatedVesicle.setNumberOfChangetypes(0);
                updatedVesicle.setNumberOfChangetracks(0);
            }
        }
    }
    else
    {
        if (currentVesicle.getVesicleId() != updatedVesicle.getVesicleId())
        {
            updatedVesicle.setNumberOfChangetypes(0);
            updatedVesicle.setNumberOfChangetracks(0);
        }

    }
    VESICLES.erase(currentVesicle);
    VESICLES.insert(updatedVesicle);
}

/*
 * Moves the Vesicle to a neighbouring microtubule
 */
void changeTrackVesicle(Vesicle v, NeighbourMicrotubule neighbourMicrotubule)
{
    Vesicle newVesicle = Vesicle(v.getVesicleId(), (v.getMicrotubuleId() + NUMBER_OF_MICROTUBULES + (1 * neighbourMicrotubule)) % NUMBER_OF_MICROTUBULES,
                                 v.getSiteId(), v.getType(), v.getNumberOfChangetypes(), v.getNumberOfChangetracks() + 1);
    updateVesicle(v, newVesicle);
    NUMBER_OF_CHANGETRACKS++;
}

/*
 * Changes the type of the Vesicle
 */
void changeTypeVesicle(Vesicle v, VesicleType vesicleType)
{
    Vesicle newVesicle = Vesicle(v.getVesicleId(), v.getMicrotubuleId(), v.getSiteId(), vesicleType, v.getNumberOfChangetypes() + 1,
                                 v.getNumberOfChangetracks());
    updateVesicle(v, newVesicle);
    NUMBER_OF_CHANGETYPES++;
}

/*
 * Determines the update and performs it
 */
void rouletteWheelUpdate(double r)
{
    double pastVesiclesPropensities = 0;
    for (set<Vesicle>::iterator it = VESICLES.begin(); it != VESICLES.end(); it++)
    {
        Vesicle v = (*it);
        Vesicle newVesicle = Vesicle(0, 0, 0, SMOOTH_ANTEROGRADE, 0, 0);
        double hoppingPropensity = getHoppingPropensity(v);
        if (hoppingPropensity > 0)
        {
            pastVesiclesPropensities += hoppingPropensity;
            if (pastVesiclesPropensities >= r)
            {
                long newVesicleId = v.getVesicleId();
                if ((v.getSiteId() == 0 && v.getDirection() == ANTEROGRADE) || (v.getSiteId() == (NUMBER_OF_SITES - 1) && v.getDirection() == RETROGRADE))
                {
                    newVesicleId = generateVesicleId();
                }
                newVesicle = Vesicle(newVesicleId, v.getMicrotubuleId(), (v.getSiteId() + NUMBER_OF_SITES + (1 * v.getDirection())) % NUMBER_OF_SITES,
                                     v.getType(), v.getNumberOfChangetypes(), v.getNumberOfChangetracks());
                updateVesicle(v, newVesicle);
                NUMBER_OF_HOPS++;
                return;
            }
            NeighbourMicrotubule neighbourMicrotubule = PREVIOUS;
            pastVesiclesPropensities += (getChangeTrackPropensity(v, neighbourMicrotubule) / SIDESTEP_FACTOR);
            if (pastVesiclesPropensities >= r)
            {
                changeTrackVesicle(v, neighbourMicrotubule);
                return;
            }
            neighbourMicrotubule = NEXT;
            pastVesiclesPropensities += (getChangeTrackPropensity(v, neighbourMicrotubule) / SIDESTEP_FACTOR);
            if (pastVesiclesPropensities >= r)
            {
                changeTrackVesicle(v, neighbourMicrotubule);
                return;
            }
            if (v.getType() == STAGGERED_ANTEROGRADE)
            {
                pastVesiclesPropensities += (StA_To_SmA_CHANGETYPE_PROPENSITY / TURNBACK_FACTOR);
                if (pastVesiclesPropensities >= r)
                {
                    changeTypeVesicle(v, SMOOTH_ANTEROGRADE);
                    return;
                }
                pastVesiclesPropensities += (StA_To_SmR_CHANGETYPE_PROPENSITY / TURNBACK_FACTOR);
                if (pastVesiclesPropensities >= r)
                {
                    changeTypeVesicle(v, SMOOTH_RETROGRADE);
                    return;
                }
            }
            else
            {
                pastVesiclesPropensities += (getChangeTypePropensity(v) / TURNBACK_FACTOR);
                if (pastVesiclesPropensities >= r)
                {
                    /* This code is to test with SmA and SmR only. For full system, remove the
                     * following 5 lines and uncomment the line with STAGGERED_ANTEROGRADE
                     */
                    /*VesicleType newType = SMOOTH_ANTEROGRADE;
                     if (v.getType() == SMOOTH_ANTEROGRADE)
                     {
                     newType = SMOOTH_RETROGRADE;
                     }
                     changeTypeVesicle(v, newType);*/
                    changeTypeVesicle(v, STAGGERED_ANTEROGRADE);
                    return;
                }
            }
        }
        else
        {
            NeighbourMicrotubule neighbourMicrotubule = PREVIOUS;
            pastVesiclesPropensities += getChangeTrackPropensity(v, neighbourMicrotubule);
            if (pastVesiclesPropensities >= r)
            {
                changeTrackVesicle(v, neighbourMicrotubule);
                return;
            }
            neighbourMicrotubule = NEXT;
            pastVesiclesPropensities += getChangeTrackPropensity(v, neighbourMicrotubule);
            if (pastVesiclesPropensities >= r)
            {
                changeTrackVesicle(v, neighbourMicrotubule);
                return;
            }

            if (v.getType() == STAGGERED_ANTEROGRADE)
            {
                pastVesiclesPropensities += StA_To_SmA_CHANGETYPE_PROPENSITY;
                if (pastVesiclesPropensities >= r)
                {
                    changeTypeVesicle(v, SMOOTH_ANTEROGRADE);
                    return;
                }
                pastVesiclesPropensities += StA_To_SmR_CHANGETYPE_PROPENSITY;
                if (pastVesiclesPropensities >= r)
                {
                    changeTypeVesicle(v, SMOOTH_RETROGRADE);
                    return;
                }
            }
            else
            {
                pastVesiclesPropensities += getChangeTypePropensity(v);
                if (pastVesiclesPropensities >= r)
                {
                    /* This code is to test with SmA and SmR only. For full system, remove the
                     * following 5 lines and uncomment the line with STAGGERED_ANTEROGRADE
                     */
                    /*VesicleType newType = SMOOTH_ANTEROGRADE;
                     if (v.getType() == SMOOTH_ANTEROGRADE)
                     {
                     newType = SMOOTH_RETROGRADE;
                     }
                     changeTypeVesicle(v, newType);*/
                    changeTypeVesicle(v, STAGGERED_ANTEROGRADE);
                    return;
                }
            }
        }
    }
}

/*
 * Returns the sum of the propensities of all vesicles.
 */
double getSumOfPropensities()
{
    double sumOfPropensities = 0;
    for (set<Vesicle>::iterator it = VESICLES.begin(); it != VESICLES.end(); it++)
    {
        sumOfPropensities += getVesiclePropensity(*it);
    }
    return (sumOfPropensities);
}

/*
 * Performs one update. Returns the time after which the update happens.
 */
double performUpdate()
{
    double sop = getSumOfPropensities();
    if (sop > RAND_MAX)
    {
        cout << "Sum of Propensities is too large" << endl;
        exit(EXIT_FAILURE);
    }
    if (sop == 0)
    {
        cout << "No more moves possible" << endl;
        exit(EXIT_FAILURE);
    }
    //long r = rand() % sop + 1; // 1 is added to avoid r taking the value 0
    double r = randomNumber();
    while (r == 0)
    {
        r = randomNumber();
        cout << "Call to rand() returns zero" << endl;
    }
    rouletteWheelUpdate(r * sop);
    double elapsedTime = (1 / sop) * log(1 / r);
    if (isinf(elapsedTime))
    {
        cout << "The value of sum of propensities is " << sop << endl;
        cout << "The random number is " << r << endl;
        exit(EXIT_FAILURE);
    }
    return (elapsedTime);
}

void writeMicrotubuleClusterDurations()
{
    ofstream durationsWriter;
    durationsWriter.open("microtubuleDurations.csv", ios::out | ios::app);
    for (int i = 0; i < NUMBER_OF_MICROTUBULES; i++)
    {
        ostringstream tmp;
        long duration = 1;
        cout << "the size of durations of microtubule " << i << " is " << MICROTUBULE_CLUSTER_DURATIONS.at(i).size() << endl;
        while ((duration < NUMBER_OF_UPDATES) && (MICROTUBULE_CLUSTER_DURATIONS.at(i).size() > 1))
        {
            long frequency = 0;
            if (MICROTUBULE_CLUSTER_DURATIONS.at(i).count(duration) > 0)
            {
                frequency = MICROTUBULE_CLUSTER_DURATIONS.at(i).at(duration);
                map<long, long>::iterator it = MICROTUBULE_CLUSTER_DURATIONS.at(i).find(duration);
                MICROTUBULE_CLUSTER_DURATIONS.at(i).erase(it);
            }
            tmp << frequency << ",";
            duration++;
        }
        durationsWriter << tmp.str() << endl;
    }
    durationsWriter.close();
}

/*
 * The values of current at a given location are written to a file.
 */
void writeCurrentToFile(vector<long> current)
{
    ofstream currentWriter;
    currentWriter.open("current.csv", ios::out | ios::app);
    ostringstream tmp;
    for (unsigned int i = 0; i < current.size(); i++)
    {
        tmp << current.at(i) << ",";
    }
    currentWriter << tmp.str() << endl;
    currentWriter.close();
}

/*
 * Output the elapsed time
 */
void writeElapsedTimeToFile()
{
    ofstream currentWriter;
    currentWriter.open("current.csv", ios::out | ios::app);
    currentWriter << fixed;
    currentWriter << TOTAL_TIME << endl;
    currentWriter.close();
}

/*
 * Output current
 */
void writeCurrent()
{
    remove("current.csv");
    writeCurrentToFile(CURRENT_AT_LOCATION_1);
    writeCurrentToFile(CURRENT_AT_LOCATION_2);
    writeCurrentToFile(CURRENT_AT_LOCATION_3);
    writeElapsedTimeToFile();
}

void testCurrentInitialization()
{
    cout << "The size of CURRENT_AT_LOCATION_1 is " << CURRENT_AT_LOCATION_1.size() << endl;
    cout << "The size of CURRENT_AT_LOCATION_2 is " << CURRENT_AT_LOCATION_2.size() << endl;
    cout << "The size of CURRENT_AT_LOCATION_3 is " << CURRENT_AT_LOCATION_3.size() << endl;
    for (unsigned int i = 0; i < CURRENT_AT_LOCATION_1.size(); i++)
    {
        cout << "The value at location " << i << " of CURRENT_AT_LOCATION_1 is " << CURRENT_AT_LOCATION_1.at(i) << endl;
        cout << "The value at location " << i << " of CURRENT_AT_LOCATION_2 is " << CURRENT_AT_LOCATION_2.at(i) << endl;
        cout << "The value at location " << i << " of CURRENT_AT_LOCATION_3 is " << CURRENT_AT_LOCATION_3.at(i) << endl;
    }
}

/*
 * Removes old files that might store current. Sets the locations that hold the current to 0
 */
void initializeCurrent()
{
    remove("current.csv");
    CURRENT_AT_LOCATION_1.assign(NUMBER_OF_TYPES, 0);
    CURRENT_AT_LOCATION_2.assign(NUMBER_OF_TYPES, 0);
    CURRENT_AT_LOCATION_3.assign(NUMBER_OF_TYPES, 0);
}

void initializeDensity()
{
    remove("density.csv");
    for (int i = 0; i < NUMBER_OF_TYPES; i++)
    {
        DENSITY.push_back(vector<long>());
    }
    for (int i = 0; i < NUMBER_OF_TYPES; i++)
    {
        DENSITY.at(i).assign(NUMBER_OF_SITES, 0);
    }
}

void initializeClusters()
{
    remove("durations.csv");
    remove("microtubuleDurations.csv");
    for (int i = 0; i < NUMBER_OF_MICROTUBULES; i++)
    {
        vector<long> microtubule;
        microtubule.assign(NUMBER_OF_SITES, -1);
        MICROTUBULE_CONTINUING_CLUSTERS.push_back(microtubule);
        map<long, long> newMap;
        MICROTUBULE_CLUSTER_DURATIONS.push_back(newMap);
    }
    CONTINUING_CLUSTERS.assign(NUMBER_OF_SITES, -1);
}

void writeDensityToFile()
{
    ofstream densityWriter;
    densityWriter.open("density.csv", ios::out | ios::app);
    for (unsigned int i = 0; i < DENSITY.size(); i++)
    {
        ostringstream tmp;
        for (unsigned int j = 0; j < DENSITY.at(i).size(); j++)
        {
            tmp << DENSITY.at(i).at(j) << ",";
        }
        densityWriter << tmp.str() << endl;
    }
    densityWriter.close();
}

void writeNumberOfUpdatesToFile()
{
    ofstream densityWriter;
    densityWriter.open("density.csv", ios::out | ios::app);
    densityWriter << THIS_UPDATE / NUMBER_OF_VESICLES << endl;
    densityWriter.close();
}

void writeDensity()
{
    writeDensityToFile();
    writeNumberOfUpdatesToFile();
}

/*
 * Write the data needed for generating kymograph
 */
void writeKymograph()
{
    int numberOfSnapshots = NUMBER_OF_SNAPSHOTS;
    double elapsedTime = 0;
    int numberOfUpdates = 0;
    while (numberOfSnapshots > 0)
    {
        while (elapsedTime < 0.2)
            //while (numberOfUpdates < (20 * UPDATES_PER_MCS))
        {
            elapsedTime += performUpdate();
            numberOfUpdates++;
        }
        writeSnapshot();
        elapsedTime = 0;
        numberOfUpdates = 0;
        numberOfSnapshots--;
    }
}

/*
 * Information on the order in which elements in VESICLES are accessed
 */
void checkIterationOrder()
{
    for (set<Vesicle>::iterator it = VESICLES.begin(); it != VESICLES.end(); it++)
    {
        Vesicle v = (*it);
        cout << "The vesicle is in microtubule number " << v.getMicrotubuleId() << " and site Number " << v.getSiteId() << endl;
    }
}

/*
 * Introduce the microtubule breaks and mitochondria in the axon
 */
/*
 * Introduce the microtubule breaks and mitochondria in the axon
 */
void makeBreaks()
{
	makeBreak(1, 574);
	makeBreak(2, 254);
	makeBreak(3, 344);
	makeBreak(4, 494);
	makeBreak(5, 494);
	makeBreak(6, 494);
	makeBreak(7, 44);
	makeBreak(8, 804);
	makeBreak(9, 624);


}

/*
 * Allow the system to settle so that measurements can be made
 */
void settleSystem()
{
    time_t beginning = time(NULL);
    for (long i = 1; i <= NUMBER_OF_UPDATES_TO_SETTLE; i++)
    {
        performUpdate();
    }
    time_t end = time(NULL);
    double settlingTime = difftime(end, beginning);
    cout << NUMBER_OF_UPDATES_TO_SETTLE << " updates were performed in " << settlingTime << " seconds" << endl;
    IS_SYSTEM_SETTLED = true;
}

/*
 *Writes the duration and frequency of clusters
 */
void writeDurations()
{
    ofstream durationsWriter;
    durationsWriter.open("durations.csv", ios::out | ios::app);
    ostringstream tmp;
    long duration = 1;
    cout << "the size of durations is " << CLUSTER_DURATIONS.size() << endl;
    while ((duration < NUMBER_OF_UPDATES) && (CLUSTER_DURATIONS.size() > 1))
    {
        long frequency = 0;
        if (CLUSTER_DURATIONS.count(duration) > 0)
        {
            frequency = CLUSTER_DURATIONS.at(duration);
            map<long, long>::iterator it = CLUSTER_DURATIONS.find(duration);
            CLUSTER_DURATIONS.erase(it);
        }
        tmp << frequency << ",";
        duration++;
    }
    durationsWriter << tmp.str() << endl;
    ostringstream tmp1;
    tmp1 << TOTAL_TIME;
    durationsWriter << tmp1.str() << endl;
    durationsWriter.close();
}

void writeMap(string filename, map<long, long> givenMap)
{
    ofstream mapWriter;
    mapWriter.open(filename.c_str(), ios::out | ios::app);
    ostringstream tmp;
    long index = 0;
    cout << "the size of the map " << filename << " is " << givenMap.size() << endl;
    while ((index < NUMBER_OF_UPDATES) && (givenMap.size() >= 1))
    {
        long frequency = 0;
        if (givenMap.count(index) > 0)
        {
            frequency = givenMap.at(index);
            map<long, long>::iterator it = givenMap.find(index);
            givenMap.erase(it);
        }
        tmp << frequency << ",";
        index++;
    }
    mapWriter << tmp.str() << endl;
    mapWriter.close();
}

void writeMapContents(string filename, map<long, long> givenMap)
{
    ofstream mapWriter;
    mapWriter.open(filename.c_str(), ios::out | ios::app);
    ostringstream tmp;
    cout << "the size of the map " << filename << " is " << givenMap.size() << endl;
    for (map<long, long>::iterator it = givenMap.begin(); it != givenMap.end(); ++it)
    {
        tmp << it->first << " => " << it->second << endl;
    }
    mapWriter << tmp.str() << endl;
    mapWriter.close();
}

void writeChangeTypeCounts()
{
    ofstream changeTypeWriter;
    changeTypeWriter.open("changetypes.csv", ios::out | ios::app);
    ostringstream tmp;
    tmp << "SmA to StA at SC is " << SmA2StAAtSC << endl;
    tmp << "SmR to StA at SC is " << SmR2StAAtSC << endl;
    tmp << "StA to SmA at SC is " << StA2SmAAtSC << endl;
    tmp << "StA to SmR at SC is " << StA2SmRAtSC << endl;
    tmp << "SmA to StA away from SC is " << SmA2StAAwayFromSC << endl;
    tmp << "SmR to StA away from SC is " << SmR2StAAwayFromSC << endl;
    tmp << "StA to SmA away from SC is " << StA2SmAAwayFromSC << endl;
    tmp << "StA to SmR away from SC is " << StA2SmRAwayFromSC << endl;
    tmp << "All SmA to StA is " << AllSmA2StA << endl;
    tmp << "All SmR to StA is " << AllSmR2StA << endl;
    tmp << "All StA to SmA is " << AllStA2SmA << endl;
    tmp << "All StA to SmR is " << AllStA2SmR << endl;
    changeTypeWriter << tmp.str() << endl;
    changeTypeWriter.close();
}

void makeMeasurements()
{
    for (THIS_UPDATE = 1; THIS_UPDATE <= NUMBER_OF_UPDATES; THIS_UPDATE++)
    {
        double currentTime = performUpdate();
        TIME_SINCE_SNAPSHOT += currentTime;
        TOTAL_TIME += currentTime;
    }
    writeCurrent();
    writeDensity();
    writeChangeTypeCounts();
    writeMapContents("durationcontents.csv", CLUSTER_DURATIONS);
    writeMapContents("changetypecontents.csv", CHANGETYPE_FREQUENCY);
    writeMapContents("changetrackcontents.csv", CHANGETRACK_FREQUENCY);
    writeDurations();
    writeMap("changetypefrequency.csv", CHANGETYPE_FREQUENCY);
    writeMap("changetrackfrequency.csv", CHANGETRACK_FREQUENCY);
    writeMicrotubuleClusterDurations();
    cout << "The elapsed time is " << TOTAL_TIME << endl;
    cout << "The number of Vesicles is " << VESICLES.size() << endl;
    cout << "The number of updates per 20 MCS is " << (NUMBER_OF_UPDATES / TOTAL_TIME) * 0.2 << endl;
    writeKymograph();
}

void initializeKymograph()
{
    for (int i = 0; i < NUMBER_OF_MICROTUBULES; i++)
    {
        ostringstream stringStream;
        stringStream << "kymo" << i << ".csv";
        string fileName;
        fileName = stringStream.str();
        remove(fileName.c_str());
    }
    for (int i = 0; i < NUMBER_OF_MICROTUBULES; i++)
    {
        ostringstream stringStream;
        stringStream << "kymoType" << i << ".csv";
        string fileName;
        fileName = stringStream.str();
        remove(fileName.c_str());
    }
}

void initializeUnblockingStatistics()
{
    remove("changetypefrequency.csv");
    remove("changetrackfrequency.csv");
}

void initializeMeasurements()
{
    initializeCurrent();
    initializeDensity();
    initializeClusters();
    initializeUnblockingStatistics();
    initializeKymograph();
}

/*
 * Display metadata on the updates that have been performed
 */
void displayMetadata()
{
    cout << "The number of hops is " << NUMBER_OF_HOPS++ << endl;
    cout << "The number of changetracks is " << NUMBER_OF_CHANGETRACKS++ << endl;
    cout << "The number of changetypes is " << NUMBER_OF_CHANGETYPES++ << endl;
    cout << "The size of CHANGETYPE_FREQUENCY is " << CHANGETYPE_FREQUENCY.size() << endl;
    cout << "The size of CHANGETRACK_FREQUENCY is " << CHANGETRACK_FREQUENCY.size() << endl;
}

void initializeAxon()
{
    makeBreaks();
    makeVesicles();
}

void printVesiclesByType()
{
    int SmACount = 0;
    int StACount = 0;
    int SmRCount = 0;
    for (set<Vesicle>::iterator it = VESICLES.begin(); it != VESICLES.end(); it++)
    {
        Vesicle v = (*it);
        if (v.getType() == SMOOTH_ANTEROGRADE)
        {
            SmACount++;
        }
        if (v.getType() == STAGGERED_ANTEROGRADE)
        {
            StACount++;
        }
        if (v.getType() == SMOOTH_RETROGRADE)
        {
            SmRCount++;
        }
    }
    cout << "The number of SmA is " << SmACount << endl;
    cout << "The number of StA is " << StACount << endl;
    cout << "The number of SmR is " << SmRCount << endl;
}

/*
 * This function initialises and settles the system and computes the current and density
 */
int main()
{
    srand(time(NULL));
    initializeAxon();
    //stuffVesicles();
    printVesiclesByType();
    settleSystem();
    initializeMeasurements();
    cout << "The size of CHANGETYPE_FREQUENCY is " << CHANGETYPE_FREQUENCY.size() << endl;
    cout << "The size of CHANGETRACK_FREQUENCY is " << CHANGETRACK_FREQUENCY.size() << endl;
    cout << "The size of VESICLES is " << VESICLES.size() << endl;
    makeMeasurements();
    displayMetadata();
    return (EXIT_SUCCESS);
}
