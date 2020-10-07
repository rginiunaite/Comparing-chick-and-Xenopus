//
// Created by giniunaite on 14/01/20.
// this code will be used for CiL and CoA using Lennard Jones potential. Langevin model will be used
//



#include "Aboria.h"
#include <Eigen/Core>
#include <algorithm> // std::unique_copy
#include <iostream>// writing on a text file
#include <fstream>
#include <math.h>
#include <assert.h>
#include <string>
#include <vector>
#include <iterator>


using namespace std;
using namespace Aboria;
using namespace Eigen; // objects VectorXd, MatrixXd



//double proportions(int n_seed) {
double proportions(int n_seed, double D, double eps_ij, double beta) {
    //VectorXi proportions(int n_seed, double beta) {
//double proportions(int n_seed, double beta) {

    bool domain_growth = false; // if false change length_x to 1100, true 300, Mayor false
    bool CiLonly = true;

    int length_x;
    if (domain_growth == false) {
        length_x = 1500;//1500; // normally 1100, Mayor, with arches 1500

    } else {
        length_x = 300;// length in x velocity of the chemoattractant matrix
    }

    int length_y = 218; // Mayor 218, ours 120;// from 218 to 130 reduction NarrowDomain 130
    double Lt_old = length_x;
    int real_length_y = 120;
    const double final_time = 1080; // 18hrs number of timesteps, 1min - 1timestep, from 6h tp 24hours. 1440 if 24hrs, 1080 - 18hrs
    double t = 0.0; // initialise time
    double dt = 0.01;//00625; // time step

    double dt_min = 0.00125;//00625; // time step
    int times = int(dt/dt_min);
    double dx = 1; // maybe 1/100
    int counter = 0; // to count simulations
    const size_t N = 5; // initial number of cells Mayor narrow domain, NarrowDomain 3
    double sigma = 2.0;
    double meanL = beta*sqrt(dt);//2;//2;//0.001;//1;//1;//0.01;//2;//0.05;//0.04; // mean movement in x velocity
    double mean = beta*sqrt(dt);//2;//2;//0.02;//2;//0.001;//1;//
    double cell_radius = 20.0;//// radius of a cell, Mayor 20.0, smallercells, ours 7.5
    double positions = cell_radius; // Mayor, only change for small cells smallercells 20.0
    const double diameter =
            2.0 * cell_radius; // diameter of a cell

    //double D = 5.0;//0.001; // diffusion coefficient for Brownian motion

    // CiL and CoA paramters, Lennard Jones
    vdouble2 sum_forces; // some of all the forces exerted on one cell, always set to zero for a new cell
    vdouble2 bound_sum_forces; // some of all the forces exerted on one cell, always set to zero for a new cell
    double npow = 4.0; // Lennard-Jones type model powers, attractive
    double mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive

    //parameters for Lennard-Jones model
    double sigma_max = diameter;//15.0; // cell diamter ( don't know why I wrote this before: maximum distance at which the cells can have an effect on each other)
    double search_param = 100.0;//100.0;//Mayor , ours 50.0 //radius to search for cell-cell interactions
    double f0 = 1.0;//1 before strength of the force
    //double eps_ij = 200.0; // depth of potential well, Lennard-Jones

    double a = 0.5; // Morse parameter
    vdouble2 force_ij; // store the value of force
    double dzdr = 0.0; // Lennard-Jones potential initialise
    vdouble2 change; // difference between the positions of two cells
    vdouble2 distchange; // for the pairwise distances calculation
    double distance = 0.0; // distance between two cells
    vdouble2 random_vector;


    double furthestCell;
    vector<double> xArray(5,0.0); // 200 when all the cells
    vector<double> xArrayOld(5,0.0);
    vector<double> yArray(5,0.0);
    vector<double> yArrayOld(5,0.0);
    vector<double> speed;
    vector<double> xpositions;
    vector<double> velocityX;
    vector<double> velocityY;
    int numberofcellsOld = N;
// cell variable
    ABORIA_VARIABLE(radius, double, "radius");
    ABORIA_VARIABLE(velocity, vdouble2, "velocity");// stores the velocity a particle moved
    ABORIA_VARIABLE(type, int, "type");// 0 if a cell is a leader, 1 if follower
    ABORIA_VARIABLE(arches, int, "arches");// 0 if has not reached arches yet, 1 if it has reached the branches
    typedef Particles<std::tuple<radius,type,arches, velocity>, 2> particle_type; // 2 stands for dimension

// will use stored value of the position of a particle
    typedef particle_type::position position;


/*
 * Domain growth section
 * */
//// domain growth parameters, domain_length
//
    double L_inf = 867.6;
    double L0 = 300;
    double t_s = 12.77;
    double alpha = 0.288;
    double k_0 = 291.2;

    VectorXd Gamma = VectorXd::Zero(length_x);
    VectorXd Gamma_old = VectorXd::Zero(length_x);

    for (int i = 0; i < length_x; i++) {
        Gamma(i) = i;//i; // change this to i *850/length_x
        Gamma_old(i) = Gamma(i);
    }


    // initialise the number of particles
    particle_type particles(10*N); // Mayor 10*N, ours N

    // initialise random number generator for particles entering the domain, appearing at the start in x and uniformly in y
    std::default_random_engine gen;
    // domain width 120
    std::uniform_real_distribution<double> uniform(cell_radius, length_y - 1 - cell_radius);// when length is 120

    //for double width (240)
    //std::uniform_real_distribution<double> uniform(cell_radius+ double(real_length_y)/2, length_y - double(real_length_y)/2 -1 - cell_radius);

    // for domain width 60
    //std::uniform_real_distribution<double> uniform(cell_radius, length_y - 1 - cell_radius);// when length is 120

    /*
* compact initialisation of particles
*/

    for (int i = 0; i < N; ++i) { //Note that here normally is N, not 1


        get<radius>(particles[i]) = cell_radius;

        // for domain width 120

            get<position>(particles[i]) = vdouble2(positions, (i + 1) * double(length_y - 1) / double(N) -
                                                                0.5 * double(length_y - 1) /
                                                                double(N)); // x=radius, uniformly in

            get<type>(particles[i]) = 0; // leaders, Mayor, comment

////// Mayor below this
            get<position>(particles[i+5]) = vdouble2(3*positions, (i + 1) * double(length_y - 1) / double(N) -
                                                                0.5 * double(length_y - 1) /
                                                                double(N)); // x=radius, uniformly in y


            get<position>(particles[i+10]) = vdouble2(5*positions, (i + 1) * double(length_y - 1) / double(N) -
                                                                0.5 * double(length_y - 1) /
                                                                double(N)); // x=radius, uniformly in y


            get<position>(particles[i+15]) = vdouble2(7*positions, (i + 1) * double(length_y - 1) / double(N) -
                                                            0.5 * double(length_y - 1) /
                                                            double(N)); // x=radius, uniformly in y

            get<position>(particles[i+20]) = vdouble2(9*positions, (i + 1) * double(length_y - 1) / double(N) -
                                                                 0.5 * double(length_y - 1) /
                                                                 double(N)); // x=radius, uniformly in y
            get<position>(particles[i+25]) = vdouble2(11*positions, (i + 1) * double(length_y - 1) / double(N) -
                                                                     0.5 * double(length_y - 1) /
                                                                     double(N)); // x=radius, uniformly in y
            get<position>(particles[i+30]) = vdouble2(13*positions, (i + 1) * double(length_y - 1) / double(N) -
                                                                     0.5 * double(length_y - 1) /
                                                                     double(N)); // x=radius, uniformly in y
            get<position>(particles[i+35]) = vdouble2(15*positions, (i + 1) * double(length_y - 1) / double(N) -
                                                                     0.5 * double(length_y - 1) /
                                                                     double(N)); // x=radius, uniformly in y
            get<position>(particles[i+40]) = vdouble2(17*positions, (i + 1) * double(length_y - 1) / double(N) -
                                                                     0.5 * double(length_y - 1) /
                                                                     double(N)); // x=radius, uniformly in y
            get<position>(particles[i+45]) = vdouble2(19*positions, (i + 1) * double(length_y - 1) / double(N) -
                                                                     0.5 * double(length_y - 1) /
                                                                     double(N)); // x=radius, uniformly in y



// for convergence test, inside the boundary

//
//            get<position>(particles[0]) = vdouble2(5*positions, (3 + 1) * double(length_y - 1) / double(N) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N)); // x=radius, uniformly in
//
//            //get<type>(particles[i]) = 0; // leaders, Mayor, comment

//////// Mayor below this
//            get<position>(particles[i]) = vdouble2(7*positions, (3 + 1) * double(length_y - 1) / double(N) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N)); // x=radius, uniformly in y


//            get<position>(particles[2]) = vdouble2(9*positions, (3 + 1) * double(length_y - 1) / double(N) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N)); // x=radius, uniformly in y
//
//
//            get<position>(particles[3]) = vdouble2(11*positions, (3 + 1) * double(length_y - 1) / double(N) -
//                                                            0.5 * double(length_y - 1) /
//                                                            double(N)); // x=radius, uniformly in y
//
//            get<position>(particles[4]) = vdouble2(13*positions, (3 + 1) * double(length_y - 1) / double(N) -
//                                            0.5 * double(length_y - 1) /double(N)); // x=radius, uniformly in y
//








////         narrow domain, NarrowDomain uncomment below
//
//        // Mayor below this
//        get<position>(particles[i+3]) = vdouble2(3*positions, (i + 1) * double(length_y - 1) / double(N) -
//                                                              0.5 * double(length_y - 1) /
//                                                              double(N)); // x=radius, uniformly in y
//
//
//        get<position>(particles[i+6]) = vdouble2(5*positions, (i + 1) * double(length_y - 1) / double(N) -
//                                                               0.5 * double(length_y - 1) /
//                                                               double(N)); // x=radius, uniformly in y
//
//
//        get<position>(particles[i+9]) = vdouble2(7*positions, (i + 1) * double(length_y - 1) / double(N) -
//                                                               0.5 * double(length_y - 1) /
//                                                               double(N)); // x=radius, uniformly in y
//
//        get<position>(particles[i+12]) = vdouble2(9*positions, (i + 1) * double(length_y - 1) / double(N) -
//                                                               0.5 * double(length_y - 1) /
//                                                               double(N)); // x=radius, uniformly in y
//        get<position>(particles[i+15]) = vdouble2(11*positions, (i + 1) * double(length_y - 1) / double(N) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N)); // x=radius, uniformly in y
//        get<position>(particles[i+18]) = vdouble2(13*positions, (i + 1) * double(length_y - 1) / double(N) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N)); // x=radius, uniformly in y
//        get<position>(particles[i+21]) = vdouble2(15*positions, (i + 1) * double(length_y - 1) / double(N) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N)); // x=radius, uniformly in y
//        get<position>(particles[i+24]) = vdouble2(17*positions, (i + 1) * double(length_y - 1) / double(N) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N)); // x=radius, uniformly in y
//        get<position>(particles[i+27]) = vdouble2(19*positions, (i + 1) * double(length_y - 1) / double(N) -
//                                                                0.5 * double(length_y - 1) /
//                                                                double(N)); // x=radius, uniformly in y
// narrow domain, NarrowDomain uncomment above





//end of Mayor



////             // for domain width 240
//            get<position>(particles[i]) = vdouble2(cell_radius,
//                                                   real_length_y / 2 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                   0.5 * double(real_length_y - 1) /
//                                                   double(N)); // x=radius, uniformly in y


//            get<position>(particles[i+5]) = vdouble2(3*cell_radius, real_length_y / 2 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                     0.5 * double(real_length_y - 1) /
//                                                     double(N)); // x=radius, uniformly in y
//
//
//            get<position>(particles[i+10]) = vdouble2(5*cell_radius,
//                                                      real_length_y / 2 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                      0.5 * double(real_length_y - 1) /
//                                                      double(N)); // x=radius, uniformly in y
//
//
//            get<position>(particles[i+15]) = vdouble2(7*cell_radius,
//                                                      real_length_y / 2 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                      0.5 * double(real_length_y - 1) /
//                                                      double(N)); // x=radius, uniformly in y
//
//            get<position>(particles[i+20]) = vdouble2(9*cell_radius,
//                                                      real_length_y / 2 + (i + 1) * double(real_length_y - 1) / double(N) -
//                                                      0.5 * double(real_length_y - 1) /
//                                                      double(N)); // x=radius, uniformly in y




//        // for domain width 60
//           get<position>(particles[i]) = vdouble2(cell_radius, (i + 1) * double(length_y - 1) / double(N) -
//                                                            0.5 * double(length_y - 1) /
//                                                            double(N)); // x=radius, uniformly in y

    }


            // check the centre of mass
            VectorXd MassCentreVector = VectorXd::Zero(7);
            double centre_of_mass =0;
            vdouble2 cellpos;

            for (int i = 0; i < particles.size(); i++){

                cellpos = get<position>(particles[i]);

                centre_of_mass += cellpos[0];

            }

            centre_of_mass = centre_of_mass/particles.size();

            MassCentreVector[0] = centre_of_mass;

            //cout << centre_of_mass << endl;
            int imass = 1; // this is to fill in the imass vector

//    // Mayor leaders and followers
//
    for (int i = 0; i < particles.size();i++ ){
        get<arches>(particles[i]) = 0;
        if (i >44){
            get<type>(particles[i]) = 0;

        }

        else{
            get<type>(particles[i]) = 1;
        }
    }




    particles.update_positions();


    // initialise neighbourhood search, note that the domain will grow in x velocity, so I initialise larger domain
    particles.init_neighbour_search(vdouble2(-20, -20), 5 * vdouble2(length_x, length_y), vbool2(false, false));


    //vtkWriteGrid("InitialCellsXenopus", t, particles.get_grid(true));

    // initialise random number generator to obtain random number between 0 and 2*pi
    std::default_random_engine gen1;
    gen1.seed( n_seed); // choose different seeds to obtain different random numbers

    //std::uniform_real_distribution<double> uniformpi(-M_PI/3.5,M_PI/3.5 +  M_PI); // 0 to up, M_PI/2 move straigth M_PI downt

    //std::normal_distribution<double> normal(M_PI/2,sigma); // normal distribution for filopodia

    // for Lennard-Jones
//    std::normal_distribution<double> normalXlead(meanL,sqrt(dt_min*1.0)); // mean specified variance 1 , actually it is standar deviation!!!!!!!!!!!
//    std::normal_distribution<double> normalX(mean,sqrt(dt_min*1.0)); // mean 0 variance 1
//    std::normal_distribution<double> normalY(0.0,sqrt(dt_min*1.0)); // mean 0 variance 1
    std::normal_distribution<double> normalXlead(meanL,1.0); // mean specified variance 1
    std::normal_distribution<double> normalX(mean,1.0); // mean 0 variance 1
    std::normal_distribution<double> normalY(0.0,1.0); // mean 0 variance 1

//initialise a random vector for each cell, entries of N(0,dtmin)

    int countfalse = 0;
    int counttrue = 0;

    int countcellsinarches = 0;


    //for each timestep
    //    while (t < final_time) {
    //while (furthestCell < 1000.0) {
      while (countcellsinarches < 41 && t < 6000.0){ //Mayor 10 if 50 cells,  NarrowDomain 6 if 30 cells
    //    while (t < 300.0){ // 300.0 for 5hours for twenty hours (1190.0) , 11hrs 610
    //   while (t < 30.0){ // 30.0 for 30min

//       while (particles.size() > 10){

        // Mayor comment this
        //      insert new cells
        if (particles.size()<50) {
        //if (counter % 100 == 0){
        bool free_position = true;
        particle_type::value_type f;

        get<position>(f) = vdouble2(cell_radius, uniform(gen)); // x=2, uniformly in y

        /*
         * loop over all neighbouring leaders within "dem_diameter" distance
         */
        for (auto tpl = euclidean_search(particles.get_query(), get<position>(f), diameter); tpl != false; ++tpl) {

            vdouble2 diffx = tpl.dx();

            if (diffx.norm() < diameter) {
                free_position = false;
                break;
            }
        }

        // all the cells are of the same type
        get<type>(f) = 1; // leaders, Mayor, comment

        if (free_position) {

            particles.push_back(f);
        }


        particles.update_positions();
        }
        // end of insert new cells
        t = t + dt;

        counter = counter + 1;
/*
      * Domain growth
      * */

        /*
    * Domain growth from here
   * */

        if (domain_growth == true) {
            //cout << "here" << endl;
            /*
             * Piecewise constant // all linear, for presentation
             * */

            Gamma(length_x - 1) = (L_inf * exp(alpha * (24.0 / final_time * t - t_s)) /
                                   (L_inf / L0 + exp(alpha * (24.0 / final_time * t - t_s)) - 1)) +
                                  k_0;

            for (int i = 0; i < length_x - 1; i++) {
                Gamma(i) = (double(i) / (length_x - 1)) * Gamma(length_x - 1);
            }



            /*
             * Domain growth
             * */


// comment this if domain does not grow
            /// update positions uniformly based on the domain growth

            vdouble2 x; // use variable x for the position of cells
            int pos;


            for (int i = 0; i < particles.size(); i++) {

                x = get<position>(particles[i]);
                // since I do not know how to do it for general case, I will do it for my specific

//            if (x[0] > Gamma(length_x - 1)) { // these are very extreme cases
//                get<position>(particles)[i] += vdouble2(Gamma(length_x - 1) - Gamma_old(length_x - 1), 0);
//            } else {

                get<position>(particles)[i] *= vdouble2((Gamma(length_x - 1)) / (Gamma_old(length_x - 1)),
                                                        1); // update position based on changes in Gamma
                //}

            }


            Gamma_old = Gamma;
// comment this if domain does not grow, up to here

//
//            // when I need to save some data of the matrix, growing
//            int counting_first = 0;
//            int counting_final = 0;
//
//            for (int a = 0; a < length_x; a++) {
//                counting_first = length_y * a;
//                counting_final = counting_first + length_y;
//                for (int k = counting_first; k < counting_final; k++) {
//                    chemo_3col(k, 0) = Gamma(a);
//                }
//            }
    // when I need to save some data of the matrix


        }
        /*
        * Domain growth to here
       * */

        /*
      * Update the position of particles
      * */


        //  create a random list of cell ids
        int check_rep = 0; // check for repetitions, 0 no rep, 1 rep



        MatrixXd positions = MatrixXd::Zero(particles.size(), 2);

        // position update
        for (int j = 0; j < particles.size(); j++) {
            vdouble2 x; // use variable x for the position of cells
            vdouble2 temppos; // temporary position of a neighbouring cell


            // old version from here
            x = get<position>(particles[j]);

            // update the position of a cell based on the deterministic and random forces exerted on it


            // deterministic force, from all the cells within the threshold distance


            sum_forces = vdouble2(0, 0);
            bound_sum_forces = vdouble2(0, 0);
            force_ij = vdouble2(0, 0);

            for (auto k = euclidean_search(particles.get_query(), x, search_param); k != false; ++k) {
                force_ij = vdouble2(0, 0);// so that if it finds itself, it would not set to previous force
                // make sure it is not the same cell
                if (get<id>(*k) != get<id>(particles[j])) {
                    change = get<position>(particles[j]) - get<position>(*k);
                    distance = change.norm();

//                    if(get<id>(*k) < N && particle_id(j) >= N){
//                        eps_ij =5;
//                    }
//                    else{
//                        eps_ij =1;
//                    }


                     npow = 4.0; // Lennard-Jones type model powers, attractive
                     mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
                    eps_ij = 4;//19.0; // 0 if no force

                    temppos = get<position>(*k);


                    if(temppos[0] >x[0] ){
                        npow = 4.0; // Lennard-Jones type model powers, attractive
                        mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
                        eps_ij =4;
                        // attraction an repulsion, sometimes comment
                        dzdr = npow * eps_ij * (2 * pow(sigma_max, mpow) / pow(distance, mpow + 1) -
                                                pow(sigma_max, npow) / pow(distance, npow + 1)); //Lennard-Jones, normally with minus but I do not put that minus when I do force_ij = f0 * dzdr * change/distance, below

                    }
                    else{
                        npow = 4.0; // Lennard-Jones type model powers, attractive
                        mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
                        eps_ij = 4;//19.0;
                        // repulsion only, if the cells only repel the cells ahead
                        dzdr = npow * eps_ij * (2 * pow(sigma_max, mpow) / pow(distance, mpow + 1));

                    }










    // IMPORTANT TO UNCOMMENT IF USUAL FORCES

//                    if (CiLonly == true) {
//                        dzdr = npow * eps_ij * (2 * pow(sigma_max, mpow)/pow(distance,mpow+1));
//                    } else {
//                        dzdr = npow * eps_ij * (2 * pow(sigma_max, mpow)/pow(distance,mpow+1) - pow(sigma_max,npow) / pow(distance,npow+1)); //Lennard-Jones
//
//                    }


                    force_ij = f0 * dzdr * change / distance;
                }

                sum_forces += force_ij;

            }

            // for the ordered sequence produce different number of random numbers, 0.00125 - 1, 0.0025 - 2, 0.005 - 4, 0.01 - 8, 0.02 - 16. In addition sqrt(1/2,4,8,..)
            //Mayor
            //random force, tryining to make the first part of this biased.
            if (get<type>(particles[j]) == 0){
//                if (t>8000){
//                    cout << "here even though should not be" << endl;
//                }
                random_vector[0] = normalXlead(gen1)+  normalXlead(gen1);// +normalXlead(gen1) +  normalXlead(gen1);//+ normalXlead(gen1) +  normalXlead(gen1) +normalXlead(gen1) +  normalXlead(gen1) +normalXlead(gen1) +  normalXlead(gen1) +normalXlead(gen1) +  normalXlead(gen1) + normalXlead(gen1) +  normalXlead(gen1) +normalXlead(gen1) +  normalXlead(gen1) ;

            }
            else{



                random_vector[0] = normalX(gen1);// + normalX(gen1) + normalX(gen1);// +normalX(gen1) + normalX(gen1)+ normalX(gen1) + normalX(gen1)+normalX(gen1) + normalX(gen1)+ normalX(gen1) + normalX(gen1)+normalX(gen1) + normalX(gen1)+ normalX(gen1) + normalX(gen1);
            }



            random_vector[1] =normalY(gen1);//+ normalY(gen1)+ normalY(gen1);//  + normalY(gen1) + normalY(gen1) + normalY(gen1) + normalY(gen1)+normalY(gen1) + normalY(gen1) + normalY(gen1) + normalY(gen1) + normalY(gen1) + normalY(gen1) + normalY(gen1) + normalY(gen1);




            x = get<position>(particles[j]) + dt * (sum_forces)+ sqrt(2.0 * D * dt) *
                                                                  (random_vector);// + bound_sum_forces); I could have random_vector/random_vector.norm()


//            if (get<type>(particles[j]) == 0){
//                cout << "interaction force " << (dt*sum_forces) << endl;
//                cout << "random " << sqrt(2.0 * dt * D)* (random_vector) << endl;
//
//            }

            // boundary forces before the position is updated

            if (x[0] < cell_radius) {
                bound_sum_forces += (cell_radius - x[0]) / x[0] * vdouble2(x[0], 0);
            }
            if (x[1] < cell_radius) {
                bound_sum_forces += (cell_radius - x[1]) / x[1] * vdouble2(0, x[1]);

            }
            if (x[1] > length_y - cell_radius - 1) {
                bound_sum_forces += (length_y - cell_radius - 1 - x[1]) / x[1] * vdouble2(0, x[1]);

            }
            if (x[0] > Gamma(length_x - 1) - cell_radius) {
                bound_sum_forces += (Gamma(length_x - 1) - cell_radius - x[0]) / x[0] * vdouble2(x[0], 0);

            }

            x = x + (bound_sum_forces); // boundary force


            //get<velocity>(particles[j]) = (sum_forces + random_vector);


            positions(j, 0) = x[0];
            positions(j, 1) = x[1];


            //if (free_position){
            //get<position>(particles[j]) = x;
            //cout << velocity.norm() << endl; // this is how much a cell moved

            //}

            // old version up to here


        }



        for (int ipos = 0; ipos < particles.size(); ipos++) {
//            if (t< 0.02){
//                cout << "previous " << get<position>(particles)[ipos] << endl;
//                }
            get<position>(particles)[ipos] = vdouble2(positions(ipos, 0), positions(ipos, 1));

        }

        // for Mayor's, delete particles greater than 850
        countcellsinarches = 0; // this has to be uncommented!!!
        for (auto p : particles) {
            if (get<position>(p)[0] > 850.0) {
                countcellsinarches = countcellsinarches + 1;
                get<arches>(p) = 1;
                 //get<alive>(p) = false;
            }
        }

//        if (countcellsinarches > 0){
//            cout << countcellsinarches << endl;
//        }


        particles.update_positions(); // not sure if needed here



        if (counter % int(60.0/dt ) == 0) { //  60/dt for every hour, 30.0/dt for every half an hour, 5 for every five minutes


        //    if (t <1078){
            //if (furthestCell > 980) { // for the one to see how long it takes till they reach the end


////            //     save at every time step
            #ifdef HAVE_VTK
                vtkWriteGrid("XenopusAheadattractBehindRepelD5eps4n4comp", t, particles.get_grid(true)); // EXPLACellsXenopusRepOnlys200D10timeBIAS18hrs
            #endif


/

        }

    }


//    /*
// * return the density of cells in domain_partition parts of the domain
// */
    const int domain_partition = int(round(Gamma(length_x -1 ) / double(55))); // number of intervals of 50 \mu m



    VectorXi proportions = VectorXi::Zero(domain_partition); // integer with number of cells in particular part

    double one_part = Gamma(length_x -1 ) / double(domain_partition);


    for (int i = 0; i < domain_partition; i++) {

        for (int j = 0; j < particles.size(); j++) {
            vdouble2 x = get<position>(particles[j]);
            if (i * one_part < x[0] && x[0] < (i + 1) * one_part) {
                proportions(i) += 1;
            }
        }

    }
  }



// parameter analysis
int main() {

    const int number_parameters = 1; // parameter range
    const int sim_num = 1;
    double eps = 4.0;

   //double vector_check_length = proportions(16); //just to know what the length is
   // cout << "ignore above" << endl;
//
    //int num_parts = vector_check_length.size(); // number of parts that I partition my domain
    //cout << "length " << num_parts << endl;
    //int num_parts = 17; // for 1800 timesteps
    //MatrixXf sum_of_all = MatrixXf::Zero(num_parts, number_parameters); // sum of the values over all simulations

////looping through D
    double D = 5.0;
//    for (int i=1; i < 6; i++){
//        if (i==1){
//            D=1.0;
//        }else{
//            D=double((i-1)*3);
//        }
//       cout << "D = " << D << endl;
//// looping through beta
    double beta = 0.0; // before was 0.2
//    for (int i=1; i < 6; i++){
//        beta = 0.01 * double(i);
//        cout << "beta = " << beta << endl;


      //    n would correspond to different seeds
        ////     parallel programming

       // VectorXd output = VectorXd::Zero(sim_num);
       double output;
 //     #pragma omp parallel for
        for (int n = 0; n < sim_num; n++) {


            //initialise the matrix to store the values
            //MatrixXi numbers = MatrixXi::Zero(num_parts, number_parameters);

            //cout << " n = " << n << endl;
            output = proportions( n, D, eps, beta);

//        // This is what I am using for MATLAB
//        ofstream output2("sepdataChickD10eps10CoACiLGROWINGDOMAIN" + to_string(n) + ".csv");
////
//        for (int i = 0; i < numbers.rows(); i++) {
//
//            for (int j = 0; j < numbers.cols(); j++) {
//
//              output2 << numbers(i, j) << ", ";
//
//                sum_of_all(i, j) += numbers(i, j);
//
//            }
//            output2 << "\n" << endl;
//        }

        }

//}// looping through D or beta






    /*
//    * will store everything in one matrix, the entries will be summed over all simulations
//    */

//
//    ofstream output3("DATAChickD10eps10CoACiLGROWINGDOMAIN.csv"); // at the end GrowingDomain
//
//
//    for (int i = 0; i < num_parts; i++) {
//
//        for (int j = 0; j < number_parameters; j++) {
//
//            output3 << sum_of_all(i, j) << ", ";
//
//        }
//        output3 << "\n" << endl;
//    }


}