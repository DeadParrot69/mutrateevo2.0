

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
std::string selstrString;//should a comma seperate vector
std::string startmutAString;//should a comma seperate vector
std::string startmutBString;//should a comma seperate vector
double mutmutrate;
double mutmuteffect;
double envmuteffect;
std::string envchangeintervalString; //should a comma seperate vector
std::string startswitchpointString; //should a comma seperate vector
std::string envtype; //should a comma seperate vector, 3 types, uniform, auto, skewed
int seed;


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
      else if (parId == "selstr") { ifs >> selstrString; }
      else if (parId == "startmutA") { ifs >> startmutAString; }
      else if (parId == "startmutB") { ifs >> startmutBString; }
      else if (parId == "mutmutrate") { ifs >> mutmutrate; }
      else if (parId == "mutmuteffect") { ifs >> mutmuteffect; }
      else if (parId == "envmuteffect") { ifs >> envmuteffect; }
      else if (parId == "envchangeinterval") { ifs >> envchangeintervalString; }
      else if (parId == "startswitchpoints") { ifs >> startswitchpointString; }
      else if (parId == "envtype") { ifs >> envtype; }

      else { std::cerr << "unknown parname in file"; exit(EXIT_FAILURE); }

    }
    else break;
  }
}

class Individual {
public:
  void set_mutA(double);
  double get_mutA();
  void set_mutB(double);
  double get_mutB();
  void set_switchpoint(double);
  double get_switchpoint();
  void set_envtrait(double);
  double get_envtrait();
private:
  double mutA;
  double mutB;
  double switchpoint;
  double envtrait;
};


void Individual::set_mutA(double MutA) {
  mutA = MutA;
};
void Individual::set_mutB(double MutB) {
  mutB = MutB;
};
void Individual::set_switchpoint(double SwitchPoint) {
  switchpoint = SwitchPoint;
};
void Individual::set_envtrait(double EnvTrait) {
  envtrait = EnvTrait;
};

double Individual::get_mutA() {
  return mutA;
};
double Individual::get_mutB() {
  return mutB;
};
double Individual::get_switchpoint() {
  return switchpoint;
};
double Individual::get_envtrait() {
  return envtrait;
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
std::vector<Individual> Initialize(int Popsize, double startmutA, double startmutB, double startswitchpoint) {
  std::vector<Individual> Population;
  for (int i = 0; i < Popsize; i++) {
    Individual NewIndividual;
    NewIndividual.set_mutA(startmutA);
    NewIndividual.set_mutB(startmutB);
    NewIndividual.set_switchpoint(startswitchpoint);
    NewIndividual.set_envtrait(0.5);
    Population.push_back(NewIndividual);
  }
  return Population;
}

//Take population object and calculate average mutation rate (+std dev) and calculate average envtrait value (+std dev) return as vector
std::vector<double> Summarise_population(std::vector<Individual> Population) {

  //Averages
  double mutrateA_sum = 0;
  double mutrateB_sum = 0;
  double envtrait_sum = 0;
  double switchpoint_sum = 0;
  for (int i = 0; i < Population.size(); ++i) {
    envtrait_sum += Population[i].get_envtrait();
    mutrateA_sum += Population[i].get_mutA();
    mutrateB_sum += Population[i].get_mutB();
    switchpoint_sum += Population[i].get_switchpoint();
  }
  double average_envtrait = envtrait_sum / static_cast <double> (Population.size());
  double average_mutrateA = mutrateA_sum / static_cast <double> (Population.size());
  double average_mutrateB = mutrateB_sum / static_cast <double> (Population.size());
  double average_switchpoint = switchpoint_sum / static_cast <double> (Population.size());

  //Standard deviations
  double mutrateAdev_sum = 0;
  double mutrateBdev_sum = 0;
  double envtraitdev_sum = 0;
  double switchpointdev_sum = 0;
  for (int i = 0; i < Population.size(); ++i) {
    envtraitdev_sum += pow((Population[i].get_envtrait() - average_envtrait), 2);
    mutrateAdev_sum += pow((Population[i].get_mutA() - average_mutrateA), 2);
    mutrateBdev_sum += pow((Population[i].get_mutB() - average_mutrateB), 2);
    switchpointdev_sum += pow((Population[i].get_switchpoint() - average_switchpoint), 2);
  }
  double std_dev_envtrait = sqrt(envtraitdev_sum / static_cast <double> (Population.size()));
  double std_dev_mutrateA = sqrt(mutrateAdev_sum / static_cast <double> (Population.size()));
  double std_dev_mutrateB = sqrt(mutrateBdev_sum / static_cast <double> (Population.size()));
  double std_dev_switchpoint = sqrt(switchpointdev_sum / static_cast <double> (Population.size()));
  std::vector<double> summary{ average_mutrateA, std_dev_mutrateA,average_mutrateB, std_dev_mutrateB, average_switchpoint, std_dev_switchpoint, average_envtrait, std_dev_envtrait };

  return summary;
}

int main(int argc, char** argv)

{
 std::mt19937_64 rng; //create rng
int seed = std::stoi(argv[argc - 1]);
  rng.seed(seed);//set seed
  
  const auto starttime = std::chrono::system_clock::now();
  std::string parameterFileName = "parameters.txt";
  readParameters(parameterFileName);

  std::cout << "maxtime: " << maxtime << '\n'
    << "seed: " << seed << '\n'
    << "saveinterval: " << saveinterval << '\n'
    << "popsize: " << popsize << '\n'
    << "couple: " << couple << '\n'
    << "selstr: " << selstrString << '\n'
    << "startmutA: " << startmutAString << '\n'
    << "startmutB: " << startmutBString << '\n'
    << "mutmutrate: " << mutmutrate << '\n'
    << "mutmuteffect: " << mutmuteffect << '\n'
    << "envmuteffect: " << envmuteffect << '\n'
    << "envchangeinterval " << envchangeintervalString << '\n'
    << "startswitchpoints " << startswitchpointString << '\n'
    << "envtype " << envtype  << std::endl;

  // turn the env change string into a vector
  std::vector<std::string> envchangevectorString = split(envchangeintervalString, ",");
  std::vector<double> envchangevector;
  for (int i = 0; i < envchangevectorString.size(); i++) {
    double nextValue = std::stod(envchangevectorString[i]);
    envchangevector.push_back(nextValue);
  }
  // turn the switch point string into a vector
  std::vector<std::string> StartSwitchpoint_vectorString = split(startswitchpointString, ",");
  std::vector<double> StartSwitchpoints;
  for (int i = 0; i < StartSwitchpoint_vectorString.size(); i++) {
    double nextValue = std::stod(StartSwitchpoint_vectorString[i]);
    StartSwitchpoints.push_back(nextValue);
  }
  // turn the mutA string into a vector
  std::vector<std::string> MutA_vectorString = split(startmutAString, ",");
  std::vector<double> MutA_Start;
  for (int i = 0; i < MutA_vectorString.size(); i++) {
    double nextValue = std::stod(MutA_vectorString[i]);
    MutA_Start.push_back(nextValue);
  }
  // turn the mutA string into a vector
  std::vector<std::string> MutB_vectorString = split(startmutBString, ",");
  std::vector<double> MutB_Start;
  for (int i = 0; i < MutB_vectorString.size(); i++) {
    double nextValue = std::stod(MutB_vectorString[i]);
    MutB_Start.push_back(nextValue);
  }
  // turn the selstr string into a vector
  std::vector<std::string> SelStr_vectorString = split(selstrString, ",");
  std::vector<double> selstrvec;
  for (int i = 0; i < SelStr_vectorString.size(); i++) {
    double nextValue = std::stod(SelStr_vectorString[i]);
    selstrvec.push_back(nextValue);
  }
  //turn the couple string into a vector
  std::vector<std::string> couplevector = split(couple, ",");
  //turn the envtype string into a vector
  std::vector<std::string> envtypevector = split(envtype, ",");




  for (int envtype_nr = 0; envtype_nr != envtypevector.size(); envtype_nr++) {//start envtype loop

    for (int selstr_nr = 0; selstr_nr != selstrvec.size(); selstr_nr++) {//start selstr loop

      for (int mutA_nr = 0; mutA_nr != MutA_Start.size(); mutA_nr++) {//start mutA loop

        for (int mutB_nr = 0; mutB_nr != MutB_Start.size(); mutB_nr++) {//start mutB loop

          for (int couple_nr = 0; couple_nr != couplevector.size(); couple_nr++) {//start couple loop

            for (int env_nr = 0; env_nr != envchangevector.size(); env_nr++) {//start envchangevector loop

              for (int switch_nr = 0; switch_nr != StartSwitchpoints.size(); switch_nr++) {//start switchpoint loop
 
                  

                  //Output file stream
                  std::ofstream ofs1("EnvType_" + envtypevector[envtype_nr] + "_Selstr_" + std::to_string(selstrvec[selstr_nr]) + "_MutA_" + std::to_string(MutA_Start[mutA_nr])
                    + "_MutB_" + std::to_string(MutB_Start[mutB_nr]) + "_Couple_" + couplevector[couple_nr] + "_EnvChange_" + std::to_string(envchangevector[env_nr]) 
                    + "_Switchpt_" + std::to_string(StartSwitchpoints[switch_nr])+"_Seed_"+ std::to_string(seed) + ".csv");

                  ofs1 << "Time" << "," << "Average_MutrateA" << "," << "Stdev_MutrateA" << "," << "Average_MutrateB" << "," << "Stdev_MutrateB" << "," << "Average_Switchpoint" << "," << "Stdev_Switchpoint" << "," << "Average_Envtrait" << "," << "Stdev_Envtrait" << "," << "Env_Value" << "," << "Seed";
                  ofs1 << '\n';


                  std::vector<Individual> Population = Initialize(popsize, MutA_Start[mutA_nr], MutB_Start[mutB_nr], StartSwitchpoints[switch_nr]);

                  //Environment distribution
                  std::uniform_real_distribution<double> uniform(0.0, 1.0);
                  std::normal_distribution<double> skewed(0.5, 0.25);



                  //Distribution mutation effect environment, mutation rates
                  std::normal_distribution<double>normalenv(0.0, envmuteffect);
                  std::normal_distribution<double>normalmut(0.0, mutmuteffect);

                  //Initialize environment - depending on env type
                  double currentenv;
                  if (envtypevector[envtype_nr] == "Skewed") { currentenv = skewed(rng); }
                  else{ currentenv = uniform(rng); }

                  int timetochange = maxtime + 1; //so it will never change in the loop...
                  if (envchangevector[env_nr] != 0) { std::geometric_distribution<int>geometric(1 / envchangevector[env_nr]); timetochange = geometric(rng); }//unless the environment is variable, in which case draw time to swtich from geometric distro
                  int switchcounter = 0;

                  //Print current paramters & time since start to console
                  const auto p1 = std::chrono::system_clock::now();
                  std::cout << "Time: " << std::chrono::duration_cast<std::chrono::seconds>(p1.time_since_epoch()).count() - std::chrono::duration_cast<std::chrono::seconds>(starttime.time_since_epoch()).count() << 
                    " Seed " << seed << ", Envchange " << envchangevector[env_nr] << ", Couple " << couplevector[couple_nr] << ", Switchpoint " << StartSwitchpoints[switch_nr] << std::endl;

                  for (int time = 0; time <= maxtime; time++) { //looping through time


                    //Check if envirionment needs changing
                    if (switchcounter == timetochange && envchangevector[env_nr] != 0) {

                      if(envtypevector[envtype_nr]=="Uniform")
                      { currentenv = uniform(rng);}

                      else if (envtypevector[envtype_nr] =="Skewed") { 
                        currentenv = skewed(rng); 
                        while ((currentenv < 0) | (currentenv>1)) {
                          currentenv = skewed(rng);
                        }
                      }
                      else if (envtypevector[envtype_nr] =="Auto") {
                        std::normal_distribution<double> skewed_auto(currentenv, 0.1);
                        currentenv = skewed_auto(rng);
                        while ((currentenv < 0) | (currentenv>1)) {
                          currentenv = skewed_auto(rng);
                        }
                      }
                      else { std::cout << "Error - env type unknown" << std::endl; }
                      


                      std::geometric_distribution<int>geometric(1 / envchangevector[env_nr]);
                      timetochange = geometric(rng);
                      switchcounter = -1;
                    }
                    //advance switchcounter for environmental change
                    switchcounter = switchcounter + 1;


                    //save the population&env
                    if (time != 1) {
                      if (time % saveinterval == 0 || time == 0) {
                        std::vector<double> summary = Summarise_population(Population);
                        ofs1 << time << "," << summary[0] << "," << summary[1] << "," << summary[2] << "," << summary[3] << "," << summary[4] << "," << summary[5] << "," << summary[6] << "," << summary[7] << "," << currentenv << "," << seed << '\n';
                      } //if remainder or time equal to zero
                    }//check time is not equal to one


                   //Calculate fitness (if environment is below or above one, then error, probability above 1 or below zero)
                    std::vector<double>fitnessvec;
                    for (int i = 0; i < Population.size(); i++) {
                      double indivfit = pow((1 - pow((currentenv - Population[i].get_envtrait()), 2.0)), selstrvec[selstr_nr]);
                      //double indivfit = 1 - abs((currentenv - Population[i].get_envtrait()));
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
                      double currentmutationrate;
                      if (fitnessvec[repvec[i]] >= (1 - Population[repvec[i]].get_switchpoint())) {
                        currentmutationrate = Population[repvec[i]].get_mutA();
                      }
                      if (fitnessvec[repvec[i]] < (1 - Population[repvec[i]].get_switchpoint())) {
                        currentmutationrate = Population[repvec[i]].get_mutB();
                      }

                      double currentmutmutationrate;
                      if (couplevector[couple_nr] == "no") { currentmutmutationrate = mutmutrate; } //If it is uncoupled use the fixed mutmutationrate
                      if (couplevector[couple_nr] == "yes") { currentmutmutationrate = currentmutationrate; } //If it is coupled use the mutationrate for that individual as its mutmutationrate

                      std::bernoulli_distribution mutmutation(currentmutmutationrate);
                      std::bernoulli_distribution envmutation(currentmutationrate);

                      Individual NewIndividual;

                      //INHERIT ENVLOCUS
                      if (envmutation(rng)) {//mutation in envlocus
                        double mutatedenv = Population[repvec[i]].get_envtrait() + normalenv(rng);
                        if (mutatedenv < 0) { mutatedenv = 0; }//Prevent env locus mutating below zero
                        if (mutatedenv > 1) { mutatedenv = 1; }//Prevent env locus mutating above one
                        NewIndividual.set_envtrait(mutatedenv);
                      }
                      else {//no mutation in envlocus
                        NewIndividual.set_envtrait(Population[repvec[i]].get_envtrait());
                      }

                      //INHERIT MUTLOCUS A
                      if (mutmutation(rng)) {//mutation in envlocus
                        double mutatedmutA = Population[repvec[i]].get_mutA() / (Population[repvec[i]].get_mutA() + (1 - Population[repvec[i]].get_mutA()) * exp(-normalmut(rng)));
                        NewIndividual.set_mutA(mutatedmutA);

                      }
                      else {//no mutation in envlocus
                        NewIndividual.set_mutA(Population[repvec[i]].get_mutA());
                      }

                      //INHERIT MUTLOCUS B
                      if (mutmutation(rng)) {//mutation in envlocus
                        double mutatedmutB = Population[repvec[i]].get_mutB() / (Population[repvec[i]].get_mutB() + (1 - Population[repvec[i]].get_mutB()) * exp(-normalmut(rng)));
                        NewIndividual.set_mutB(mutatedmutB);

                      }
                      else {//no mutation in envlocus
                        NewIndividual.set_mutB(Population[repvec[i]].get_mutB());
                      }

                      //INHERIT Switchpoint
                      if (mutmutation(rng)) {//mutation in envlocus
                        double mutatedswitchpoint = Population[repvec[i]].get_switchpoint() / (Population[repvec[i]].get_switchpoint() + (1 - Population[repvec[i]].get_switchpoint()) * exp(-normalmut(rng)));
                        NewIndividual.set_switchpoint(mutatedswitchpoint);

                      }
                      else {//no mutation in envlocus
                        NewIndividual.set_switchpoint(Population[repvec[i]].get_switchpoint());
                      }



                      // NewIndividual.set_switchpoint(Population[repvec[i]].get_switchpoint());

                       //put new individual in population
                      NewPop.push_back(NewIndividual);

                    }//end inheritance forloop

                    Population = NewPop;//update population


                  }//end time loop
              }//end switch loop
            }//end envchange loop
          }//end couple loop
        }//end mutB loop
      }//end mutA loop
    }//end selstr loop
  }//end envtype loop

}