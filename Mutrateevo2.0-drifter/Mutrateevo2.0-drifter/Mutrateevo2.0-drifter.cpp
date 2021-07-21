

#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <string>
#include <iomanip>
#include <utility>
#include <chrono>

int maxtime;
int saveinterval;
int popsize;
std::string couple; //value yes or no, seperated by comma
std::string startmutAString;//should a comma seperate vector
double mutmutrate;
double mutmuteffect;




//Reading in parameters from file
void readParameters(const std::string& parFileName)


{
  std::ifstream ifs(parFileName.c_str());
  if (!ifs.is_open()) { std::cerr << "Unable to open parfile " << parFileName << '\n'; exit(EXIT_FAILURE); }
  for (;;) {
    std::string parId;
    ifs >> parId;
    if (ifs.good()) {
      if (parId == "maxtime") { ifs >> maxtime; }
      else if (parId == "saveinterval") { ifs >> saveinterval; }
      else if (parId == "popsize") { ifs >> popsize; }
      else if (parId == "couple") { ifs >> couple; }
      else if (parId == "startmutA") { ifs >> startmutAString; }
      else if (parId == "mutmutrate") { ifs >> mutmutrate; }
      else if (parId == "mutmuteffect") { ifs >> mutmuteffect; }
   

      else { std::cerr << "unknown parname in file"; exit(EXIT_FAILURE); }

    }
    else break;
  }
}

class Individual {
public:
  void set_mutA(double);
  double get_mutA();
private:
  double mutA;
};


void Individual::set_mutA(double MutA) {
  mutA = MutA;
};


double Individual::get_mutA() {
  return mutA;
};



//splits a string "str" each time the token appears
std::vector<std::string> split(std::string str, std::string token) {
  std::vector<std::string> result;
  while (str.size()) {
    int index = str.find(token);
    if (index != std::string::npos) {
      result.push_back(str.substr(0, index));
      str = str.substr(index + token.size());
      if (str.size() == 0)result.push_back(str);
    }
    else {
      result.push_back(str);
      str = "";
    }
  }
  return result;
}

//Initializing the population
std::vector<Individual> Initialize(int Popsize, double startmutA) {
  std::vector<Individual> Population;
  for (int i = 0; i < Popsize; i++) {
    Individual NewIndividual;
    NewIndividual.set_mutA(startmutA);
    Population.push_back(NewIndividual);
  }
  return Population;
}

//Take population object and calculate average mutation rate (+std dev) and calculate average envtrait value (+std dev) return as vector
std::vector<double> Summarise_population(std::vector<Individual> Population) {

  //Averages
  double mutrateA_sum = 0;
  for (int i = 0; i < Population.size(); ++i) {
    mutrateA_sum += Population[i].get_mutA();

  }

  double average_mutrateA = mutrateA_sum / static_cast <double> (Population.size());

  //Standard deviations
  double mutrateAdev_sum = 0;
  for (int i = 0; i < Population.size(); ++i) {
    mutrateAdev_sum += pow((Population[i].get_mutA() - average_mutrateA), 2);
  }
  double std_dev_mutrateA = sqrt(mutrateAdev_sum / static_cast <double> (Population.size()));

  std::vector<double> summary{ average_mutrateA, std_dev_mutrateA};

  return summary;
}

int main(int argc, char** argv)

{
  std::mt19937_64 rng; //create rng
  int seed = std::stoi(argv[argc - 1]);
  rng.seed(seed);//set seed

  const auto starttime = std::chrono::system_clock::now();
  std::string parameterFileName = "parameters_drifter.txt";
  readParameters(parameterFileName);

  std::cout << "maxtime: " << maxtime << '\n'
    << "seed: " << seed << '\n'
    << "saveinterval: " << saveinterval << '\n'
    << "popsize: " << popsize << '\n'
    << "couple: " << couple << '\n'
    << "startmutA: " << startmutAString << '\n'
    << "mutmutrate: " << mutmutrate << '\n'
    << "mutmuteffect: " << mutmuteffect << std::endl;


  // turn the mutA string into a vector
  std::vector<std::string> MutA_vectorString = split(startmutAString, ",");
  std::vector<double> MutA_Start;
  for (int i = 0; i < MutA_vectorString.size(); i++) {
    double nextValue = std::stod(MutA_vectorString[i]);
    MutA_Start.push_back(nextValue);
  }


  //turn the couple string into a vector
  std::vector<std::string> couplevector = split(couple, ",");



      for (int mutA_nr = 0; mutA_nr != MutA_Start.size(); mutA_nr++) {//start mutA loop

        for (int couple_nr = 0; couple_nr != couplevector.size(); couple_nr++) {//start couple loop






                //Output file stream
            std::ofstream ofs1("MutA_" + std::to_string(MutA_Start[mutA_nr])
              + "_Couple_" + couplevector[couple_nr] 
              + "_Seed_" + std::to_string(seed) + ".csv");

            ofs1 << "Time" << "," << "Average_MutrateA" << "," << "Stdev_MutrateA" << ","  << "Seed";
            ofs1 << '\n';


            std::vector<Individual> Population = Initialize(popsize, MutA_Start[mutA_nr]);

         
            //Distribution mutation effect environment, mutation rates
            std::normal_distribution<double>normalmut(0.0, mutmuteffect);

           

            //Print current paramters & time since start to console
            const auto p1 = std::chrono::system_clock::now();
            std::cout << "Time: " << std::chrono::duration_cast<std::chrono::seconds>(p1.time_since_epoch()).count() - std::chrono::duration_cast<std::chrono::seconds>(starttime.time_since_epoch()).count() <<
              " Seed " << seed  << ", Couple " << couplevector[couple_nr] << std::endl;

            for (int time = 0; time <= maxtime; time++) { //looping through time





              //save the population&env
              if (time != 1) {
                if (time % saveinterval == 0 || time == 0) {
                  std::vector<double> summary = Summarise_population(Population);
                  ofs1 << time << "," << summary[0] << "," << summary[1] << ","  << seed << '\n';
                } //if remainder or time equal to zero
              }//check time is not equal to one


             //Calculate fitness (if environment is below or above one, then error, probability above 1 or below zero)
              std::vector<double>fitnessvec;
              for (int i = 0; i < Population.size(); i++) {
                double indivfit = 0.5;
                fitnessvec.push_back(indivfit);
              }

              //Draw who will reproduce
              std::vector<int> repvec; //vector that contains positions of reproducing individual
              std::discrete_distribution<int> weightedLottery(fitnessvec.begin(), fitnessvec.end());
              for (int i = 0; i < popsize; ++i) {
                repvec.push_back(weightedLottery(rng));
              }

              std::vector<Individual> NewPop;//make an empty vector of individuals, this will be the new population

              //Add new individuals, with mutatation
              for (int i = 0; i < repvec.size(); i++) {
                double currentmutationrate = Population[repvec[i]].get_mutA();


                double currentmutmutationrate;
                if (couplevector[couple_nr] == "no") { currentmutmutationrate = mutmutrate; } //If it is uncoupled use the fixed mutmutationrate
                if (couplevector[couple_nr] == "yes") { currentmutmutationrate = currentmutationrate; } //If it is coupled use the mutationrate for that individual as its mutmutationrate

                std::bernoulli_distribution mutmutation(currentmutmutationrate);
    

                Individual NewIndividual;

                
                //INHERIT MUTLOCUS A
                if (mutmutation(rng)) {//mutation in envlocus
                  double mutatedmutA = Population[repvec[i]].get_mutA() / (Population[repvec[i]].get_mutA() + (1 - Population[repvec[i]].get_mutA()) * exp(-normalmut(rng)));
                  NewIndividual.set_mutA(mutatedmutA);

                }
                else {//no mutation in envlocus
                  NewIndividual.set_mutA(Population[repvec[i]].get_mutA());
                }

                //put new individual in population
                NewPop.push_back(NewIndividual);

              }//end inheritance forloop

              Population = NewPop;//update population


            }//end time loop
        }//end couple loop
      }//end mutA loop


}