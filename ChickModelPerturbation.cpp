//
// Created by rasa on 17/06/2020.
// Model for the chick host cells in perturbation experiments
//




#include "Aboria.h"
#include <Eigen/Core>
#include <algorithm> // std::unique_copy
#include <iostream>// writing on a text file
#include <fstream>
#include <math.h>
#include <assert.h>


using namespace std;
using namespace Aboria;
using namespace Eigen; // objects VectorXd, MatrixXd


VectorXi proportions(int n_seed, double D, double eps_ij, double beta) {
//double proportions(int n_seed, double D, double eps_ij, double beta) {
//VectorXi proportions(int n_seed, double beta) {
//double proportions(int n_seed, double beta) {

    bool domain_growth = true; // if false change length_x to 1100, true 300, Mayor false
    bool CiLonly = false;

    int length_x;
    if (domain_growth == false) {
        length_x = 1100;//1500; // normally 1100, Mayor, with arches 1500

    } else {
        length_x = 300;// length in x velocity of the chemoattractant matrix
    }

    int length_y = 120; // Mayor 218, ours 120;// from 218 to 130 reduction NarrowDomain 130
    double Lt_old = length_x;
    int real_length_y = 120;
    const double final_time = 1080; // 18hrs number of timesteps, 1min - 1timestep, from 6h tp 24hours. 1440 if 24hrs, 1080 - 18hrs
    double t = 0.0; // initialise time
    double dt = 0.005; // time step
    double abl_time = 300*int( 1/dt);
    int insert_time =600 *int( 1/dt);// time when new cells are inserted from donor to host, times 1/dt the hours because it checks with counter, 240 *int( 1/dt); for front-to-back and back-to-front experiments. 600 *int( 1/dt) for Removal of cells at 200-300 \mu m
    double dx = 1; // maybe 1/100
    int counter = 0; // to count simulations
    const size_t N = 5; // initial number of cells Mayor narrow domain, NarrowDomain 3
    double sigma = 2.0;
    double meanL = beta*sqrt(dt);//0.007;//beta;//0.02 // mean movement in x velocity for leaders
    double mean = 0;//beta*sqrt(dt);//beta; //0.02; // mean for all other cells
    double cell_radius = 7.5;//// radius of a cell, Mayor 20.0, smallercells, ours 7.5
    double positions = cell_radius; // Mayor, only change for small cells smallercells 20.0
    const double diameter =
            2.0 * cell_radius; // diameter of a cell
    double eps = 1; // for phenotypic switching, the distance has to be that much higher

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
    //double eps_ij = 1.0; // depth of potential well, Lennard-Jones

    double a = 0.5; // Morse parameter
    vdouble2 force_ij; // store the value of force
    double dzdr = 0.0; // Lennard-Jones potential initialise
    vdouble2 change; // difference between the positions of two cells
    vdouble2 distchange; // for the pairwise distances calculation
    double distance = 0.0; // distance between two cells
    vdouble2 random_vector;


    /*
     * FRONT TO BACK
    */
//    // positions of inserted cells
//    // ahead attract, behind repel
//    vdouble2 cell1pos = vdouble2(217.706 - 217.706 - cell_radius,26.6998);
//    vdouble2 cell2pos = vdouble2(188.849- 217.706 - cell_radius,27.2011);
//    vdouble2 cell3pos = vdouble2(196.541- 217.706 - cell_radius,98.4631);
//    vdouble2 cell4pos = vdouble2(178.204- 217.706 - cell_radius,90.29);
//    vdouble2 cell5pos = vdouble2(172.781- 217.706 - cell_radius,7.5);


//    // leaders and followers, ahead attract, behind do nothing
//    vdouble2 cell1pos = vdouble2(147.727-147.727-cell_radius,72.9717);
//    vdouble2 cell2pos = vdouble2(128.703-147.727-cell_radius,76.2494);
//    vdouble2 cell3pos = vdouble2(111.338-147.727-cell_radius,67.3626);
//    vdouble2 cell4pos = vdouble2(107.058-147.727-cell_radius,85.7495);
//    vdouble2 cell5pos = vdouble2(100.038-147.727-cell_radius,50.7496);

    // for the repulsion only model
//    vdouble2 cell1pos = vdouble2(298.219 - 298.219 - cell_radius,24.8353);
//    vdouble2 cell2pos = vdouble2(274.003- 298.219 - cell_radius,66.1192);
//    vdouble2 cell3pos = vdouble2(297.314- 298.219 - cell_radius,83.0196);
//    vdouble2 cell4pos = vdouble2( 259.136- 298.219 - cell_radius,10.1794);
//    vdouble2 cell5pos = vdouble2(199.867- 298.219 - cell_radius,28.1161);


    /*
     * BACK TO FRONT
    */
    // positions of inserted cells furthest position + 50
    // ahead attract, behind repel
//    vdouble2 cell1pos = vdouble2(14.0318,40.0887);
//    vdouble2 cell2pos = vdouble2(11.2875,11.375);
//    vdouble2 cell3pos = vdouble2(9.51958,68.74);
//    vdouble2 cell4pos = vdouble2(7.6427,29.1878);
//    vdouble2 cell5pos = vdouble2( 7.69368,52.114);


//    // leaders and followers, ahead attract, behind do nothing, PositionsLead
    vdouble2 cell1pos = vdouble2(7.5,73.4823); // adding 120, so that it would similar to the ahead attract, behind repel, could think more where they should be placed
    vdouble2 cell2pos = vdouble2(7.60203,7.5);
    vdouble2 cell3pos = vdouble2(7.5,29.0781);
    vdouble2 cell4pos = vdouble2(7.50454,97.8988);
    vdouble2 cell5pos = vdouble2(7.5 ,83.2659);


    //////    //Repulsion only, we choose from PositionsReponlyeps19D12at4h.csv at Perturbations data positions with x less than 11, max at 275.645 +100
//    vdouble2 cell1pos = vdouble2(10.314,58.9851);
//    vdouble2 cell2pos = vdouble2(9.44441,30.9005);
//    vdouble2 cell3pos = vdouble2(10.0516,7.63077);
//    vdouble2 cell4pos = vdouble2(17.4627,76.3048);
//    vdouble2 cell5pos = vdouble2(7.5,90.7541);

    // old version
    // positions of inserted cells furthest position + 50
    // ahead attract, behind repel
//    vdouble2 cell1pos = vdouble2(14.0318+ 217.706+ 50,40.0887);
//    vdouble2 cell2pos = vdouble2(11.2875+ 217.706+ 50,11.375);
//    vdouble2 cell3pos = vdouble2(9.51958+ 217.706+ 50,68.74);
//    vdouble2 cell4pos = vdouble2(7.6427+ 217.706+ 50,29.1878);
//    vdouble2 cell5pos = vdouble2( 7.69368+ 217.706+ 50,52.114);

//    // leaders and followers, ahead attract, behind do nothing, PositionsLead
//    vdouble2 cell1pos = vdouble2(7.5 + 147.727 +120,73.4823); // adding 120, so that it would similar to the ahead attract, behind repel, could think more where they should be placed
//    vdouble2 cell2pos = vdouble2(7.60203+ 147.727 +120,7.5);
//    vdouble2 cell3pos = vdouble2(7.5 + 147.727 +120,29.0781);
//    vdouble2 cell4pos = vdouble2(7.50454 + 147.727 +120,97.8988);
//    vdouble2 cell5pos = vdouble2(7.5 + 147.727 +120,83.2659);


//////    //Repulsion only, we choose from PositionsReponlyeps19D12at4h.csv at Perturbations data positions with x less than 11, max at 275.645 +100
//    vdouble2 cell1pos = vdouble2(10.314 + 298.219,58.9851);
//    vdouble2 cell2pos = vdouble2(9.44441+ 298.219,30.9005);
//    vdouble2 cell3pos = vdouble2(10.0516+ 298.219,7.63077);
//    vdouble2 cell4pos = vdouble2(17.4627+ 298.219,76.3048);
//    vdouble2 cell5pos = vdouble2(7.5+ 298.219,90.7541);




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

//// when I need to save some data of the matrix
    MatrixXd chemo_3col(length_x * length_y, 4), chemo_3col_ind(length_x * length_y,
                                                                2); // need for because that is how paraview accepts data, third dimension is just zeros
    MatrixXd chemo = MatrixXd::Zero(length_x, length_y);

    for (int i = 0; i < length_x; i++) {
        for (int j = 0; j < length_y; j++) {
            chemo(i, j) = 1;//cos(
            // M_PI * i / space_grid_controller);//1;//C0 - 0.5 * cos( M_PI * i/space_grid_controller * n);
        }
    }


    // x, y coord, 1st and 2nd columns respectively
    int k = 0;
    // it has to be 3D for paraview
    while (k < length_x * length_y) {
        for (int i = 0; i < length_x; i++) {
            for (int j = 0; j < length_y; j++) {
                chemo_3col_ind(k, 0) = i;
                chemo_3col_ind(k, 1) = j;
                chemo_3col(k, 2) = 0;
                k += 1;
            }
        }
    }

    // save the x coordinates, scaling only based on the grid
    for (int i = 0; i < length_x * length_y; i++) {
//        chemo_3col(i, 0) = Gamma_initial * chemo_3col_ind(i, 0) / double(space_grid_controller);
        chemo_3col(i, 0) = chemo_3col_ind(i, 0);// / double(space_grid_controller);

    }



    // u column
    for (int i = 0; i < length_x * length_y; i++) {
        chemo_3col(i, 3) = chemo(chemo_3col_ind(i, 0), chemo_3col_ind(i, 1));
    }


    // y coordinates, 1D so nothing changes
    for (int i = 0; i < length_x * length_y; i++) {
        // chemo_3col(i, 1) = y_init * chemo_3col_ind(i, 1) / double(space_grid_controller);
        chemo_3col(i, 1) = chemo_3col_ind(i, 1);// / double(space_grid_controller);
    }


    // to save the matrix
            int counting_first = 0;
            int counting_final = 0;

            for (int a = 0; a < length_x; a++) {
                counting_first = length_y * a;
                counting_final = counting_first + length_y;
                for (int k = counting_first; k < counting_final; k++) {
                    chemo_3col(k, 0) = Gamma(a);
                }
            }
// when I need to save some data of the matrix


    // initialise the number of particles
    particle_type particles(N); // Mayor 10*N, ours N

    // initialise random number generator for particles entering the domain, appearing at the start in x and uniformly in y
    std::default_random_engine gen;
    // domain width 120
    std::uniform_real_distribution<double> uniform(cell_radius, length_y - 1 - cell_radius);// when length is 120

/*
* compact initialisation of particles
*/

    for (int i = 0; i < N; ++i) { // initial 1 cell


        get<radius>(particles[i]) = cell_radius;

        // for domain width 120

        // one cell in the middle
        get<position>(particles[i]) = vdouble2(positions, (i + 1) * double(length_y - 1) / double(N) -
                                                          0.5 * double(length_y - 1) /
                                                          double(N)); // x=radius, uniformly in

        get<type>(particles[i]) = 0; // leaders, Mayor, comment

    }





    particles.update_positions();


    // initialise neighbourhood search, note that the domain will grow in x velocity, so I initialise larger domain
    particles.init_neighbour_search(vdouble2(-110, -200),  vdouble2(4 *length_x, 1.5*length_y), vbool2(false, false));


    //vtkWriteGrid("InitialCells", t, particles.get_grid(true));

    // initialise random number generator to obtain random number between 0 and 2*pi
    std::default_random_engine gen1;
    gen1.seed( n_seed); // choose different seeds to obtain different random numbers

    //std::uniform_real_distribution<double> uniformpi(-M_PI/3.5,M_PI/3.5 +  M_PI); // 0 to up, M_PI/2 move straigth M_PI downt

    //std::normal_distribution<double> normal(M_PI/2,sigma); // normal distribution for filopodia

    // for Lennard-Jones
    std::normal_distribution<double> normalXlead(meanL,1.0); // mean specified variance 1
    std::normal_distribution<double> normalX(mean,1.0); // mean 0 variance 1
    std::normal_distribution<double> normalY(0.0,1.0); // mean 0 variance 1

    int countfalse = 0;
    int counttrue = 0;

    int countcellsinarches = 0;



    //for each timestep
    while (t < final_time) {//  && countcellsinarches < 5) { // the second part is to check time to invasion: for whether invasion does not happen too early!!!
        // while (t < 480.0) {//  && countcellsinarches < 5) { // this is for the perturbation experiment, cut from the front or back from the donor and paste it back.

        //while (furthestCell < 1000.0) {
        //     while (countcellsinarches < 41 && t < 3000.0){ //Mayor 10 if 50 cells,  NarrowDomain 6 if 30 cells
        //  while (t < 1190.0){ // for twenty hours
//       while (particles.size() > 10){

        // Mayor comment this
        //      insert new cells
        // BACK TO FRONT
        //if (particles.size() < 50) {
        //if (counter < insert_time) { // this is for the perturbation experiments, FRONT TO BACK because the influx is ceased once new cells are inserted,
            if (counter < abl_time){// Ablation
//              if (counter % 800 == 0){
            //if (counter % 50 == 0){
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



        // insert some cells at the neural tube, just before the domain
        // very inefficient
//        if (counter == insert_time){//48000){
//            // BACK TO FRONT uncomment, positions based on domain length and the position of the most forward cell
//            vector<double> vec {0,0,0,0,0};
//            for (int i = 0; i<N; i++){
//                vec[0]=  get<position>(particles[i])[0] ;
//            }
//            double max = *max_element(vec.begin(), vec.end());
//
//
//            cell1pos[0] = cell1pos[0] + (Gamma(length_x - 1)+ max)/2; // half way between furthest cell and domain length plus the position of those cells
//            cell2pos[0] = cell2pos[0] + (Gamma(length_x - 1)+ max)/2;
//            cell3pos[0] = cell3pos[0] + (Gamma(length_x - 1)+ max)/2;
//            cell4pos[0] = cell4pos[0] + (Gamma(length_x - 1)+ max)/2;
//            cell5pos[0] = cell5pos[0] + (Gamma(length_x - 1)+ max)/2;
//
//
//            particle_type::value_type f;
//
//            get<position>(f) = cell1pos; // x=2, uniformly in y
//
//            // all the cells are of the same type
//            get<type>(f) = 1; // leaders, Mayor, comment
//
//            particles.push_back(f);
//
//            particle_type::value_type f1;
//
//            get<position>(f1) = cell2pos; // x=2, uniformly in y
//
//            // all the cells are of the same type
//            get<type>(f1) = 1; // leaders, Mayor, comment
//
//            particles.push_back(f1);
//
//            particle_type::value_type f2;
//
//            get<position>(f2) = cell3pos; // x=2, uniformly in y
//
//            // all the cells are of the same type
//            get<type>(f2) = 1; // leaders, Mayor, comment
//
//            particles.push_back(f2);
//
//            particle_type::value_type f3;
//
//            get<position>(f3) = cell4pos; // x=2, uniformly in y
//
//            // all the cells are of the same type
//            get<type>(f3) = 1; // leaders, Mayor, comment
//
//            particles.push_back(f3);
//
//            particle_type::value_type f4;
//
//            get<position>(f4) = cell5pos; // x=2, uniformly in y
//
//            // all the cells are of the same type
//            get<type>(f4) = 1; // leaders, Mayor, comment
//
//            particles.push_back(f4);
//
//
//            particles.update_positions();
//
////            for (auto p : particles) {
////                cout << get<position>(p) << endl;
////            }
//
//        }
//

//
//        // REMOVAL OF CELLS IN THE MIDDLE OF THE STREAM, AT 10H, REMOVE FROM 200-300\MU M
//        if (counter == insert_time) {//48000){
//            for (auto p : particles) {
//                if (get<position>(p)[0] > 200.0 && get<position>(p)[0] < 300.0) {
//                    get<alive>(p) = false;
//                }
//            }
//        }




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
                force_ij = vdouble2(0, 0); // if it finds itself, the force has to be 0
                // make sure it is not the same cell
                if (get<id>(*k) != get<id>(particles[j])) {
                    change = get<position>(particles[j]) - get<position>(*k);
                    distance = change.norm();
                    // if leaders attract followers more than followers leaders and leaders leaders

//                    // forces for n=2, eps =51,102,153,204,256
//                    eps_ij = 51.0;//153;//51.0;//56.0;
//                     npow = 2.0; // Lennard-Jones type model powers, attractive
//                     mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
//
//                    if(get<type>(*k) == 0 && j >= get<type>(particles[j]) == 1){//LEADERS attract followers
//                        eps_ij =51.0;//153.0;//94.0;
//
//                        npow = 2.0; // Lennard-Jones type model powers, attractive
//                        mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
//
//                    }
//                    else{
//                        eps_ij = 51.0;//153;//51.0;//56.0;
//
//                        npow = 2.0; // Lennard-Jones type model powers, attractive
//                        mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
//
//                    }


//// based on position if ahead, larger force than if behind

//                    npow = 4.0; // Lennard-Jones type model powers, attractive
//                    mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
//                    eps_ij = 1.5;//19.0; // 0 if no force
//
//                    temppos = get<position>(*k);
//
//
//                    if(temppos[0] >x[0] ){
//                        npow = 4.0; // Lennard-Jones type model powers, attractive
//                        mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
//                        eps_ij =1.5;
//                        // attraction an repulsion, sometimes comment
//                        dzdr = npow * eps_ij * (2 * pow(sigma_max, mpow) / pow(distance, mpow + 1) -
//                                                pow(sigma_max, npow) / pow(distance, npow + 1)); //Lennard-Jones, normally with minus but I do not put that minus when I do force_ij = f0 * dzdr * change/distance, below
//
//                    }
//                    else{
//                        npow = 4.0; // Lennard-Jones type model powers, attractive
//                        mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
//                        eps_ij = 1.5;//19.0;
//                        // repulsion only, if the cells only repel the cells ahead
//                        dzdr = npow * eps_ij * (2 * pow(sigma_max, mpow) / pow(distance, mpow + 1));
//
//                    }

//// based on position if ahead, no force ad force, leaders and followers

                    npow = 4.0; // Lennard-Jones type model powers, attractive
                    mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
                    eps_ij = 0.0;//19.0; // 0 if no force

                    temppos = get<position>(*k);


                    if(temppos[0] >x[0] ){
                        npow = 2.0; // Lennard-Jones type model powers, attractive
                        mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
                        eps_ij =256.0;
                        // attraction an repulsion, sometimes comment
                        dzdr = npow * eps_ij * (2 * pow(sigma_max, mpow) / pow(distance, mpow + 1) -
                                                pow(sigma_max, npow) / pow(distance, npow + 1)); //Lennard-Jones, normally with minus but I do not put that minus when I do force_ij = f0 * dzdr * change/distance, below

                    }
                    else{
                        npow = 4.0; // Lennard-Jones type model powers, attractive
                        mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
                        eps_ij = 0.0;//19.0;
                        // repulsion only, if the cells only repel the cells ahead
                        dzdr = npow * eps_ij * (2 * pow(sigma_max, mpow) / pow(distance, mpow + 1));

                    }


//                     // if gradually decreasing force
////
//                     npow = 2.0; // Lennard-Jones type model powers, attractive
//                     mpow = 2.0*npow; //  Lennard-Jones type model powers, repulsive
//                     if (get<id>(*k) < 5) {
//
//                         eps_ij = 102;//94.0;
//                     }
//                     if (get<id>(*k) < 10 && get<id>(*k) >= 5) {
//
//                         eps_ij = 51;//75.0;
//                     }
//                     if (get<id>(*k) < 15 && get<id>(*k) >= 10) {
//
//                         eps_ij = 51;//56.0;
//                     }
//                     if (get<id>(*k) < 20 && get<id>(*k) >= 15) {
//
//                         eps_ij = 51;//38.0;
//                     }
//                    if (get<id>(*k) >= 20 ){
//
//                        eps_ij = 51.0;
//                    }

//                     if (CiLonly == true) {
//                         dzdr = npow * eps_ij * (2 * pow(sigma_max, mpow) / pow(distance, mpow + 1));
//                     } else {
//                         dzdr = npow * eps_ij * (2 * pow(sigma_max, mpow) / pow(distance, mpow + 1) -
//                                                 pow(sigma_max, npow) / pow(distance, npow + 1)); //Lennard-Jones, normally with minus but I do not put that minus when I do force_ij = f0 * dzdr * change/distance, below
//
//                     }


                    force_ij = f0 * dzdr * change / distance;
                    //cout << "force_ij " << force_ij << endl;
                }

                sum_forces += force_ij;

            }


            //Mayor
            //random force, tryining to make the first part of this biased.
            if (get<type>(particles[j]) == 0) {
//                if (t>8000){
//                    cout << "here even though should not be" << endl;
//                }
                random_vector[0] = normalXlead(gen1); // 0.05;//
                // random_vector[0] = 0.05;//


            } else {
                random_vector[0] =normalX(gen1); // 0.0;//
                // random_vector[0] =0;//normalX(gen1); // 0.0;//
            }


            //random_vector[0] = normalX(gen1);
            random_vector[1] = normalY(gen1);



            x = get<position>(particles[j]) + dt * (sum_forces) + sqrt(2.0 * dt * D) *
                                                                  (random_vector);// + bound_sum_forces); I could have random_vector/random_vector.norm()


//            if (get<type>(particles[j]) == 0){
//                cout << "interaction force " << (dt*sum_forces) << endl;
//                cout << "random " << sqrt(2.0 * dt * D)* (random_vector) << endl;
//
//            }
//
//             boundary forces before the position is updated

            if (x[0] < cell_radius) {
                bound_sum_forces += (cell_radius- x[0]) / x[0] * vdouble2(x[0], 0);
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

//           //  FORCES change for Front to BACK experiments
//           //  allow till -100
//            if (counter < insert_time){ // this is before the insertion, domain as usual
//                if (x[0] < cell_radius) {
//                    bound_sum_forces += (cell_radius- x[0]) / x[0] * vdouble2(x[0], 0);
//                }
//                if (x[1] < cell_radius) {
//                    bound_sum_forces += (cell_radius - x[1]) / x[1] * vdouble2(0, x[1]);
//
//                }
//                if (x[1] > length_y - cell_radius - 1) {
//                    bound_sum_forces += (length_y - cell_radius - 1 - x[1]) / x[1] * vdouble2(0, x[1]);
//
//                }
//                if (x[0] > Gamma(length_x - 1) - cell_radius) {
//                    bound_sum_forces += (Gamma(length_x - 1) - cell_radius - x[0]) / x[0] * vdouble2(x[0], 0);
//
//                }
//
//                x = x + (bound_sum_forces); // boundary force
//            }
//            else{
//                if (x[0] < cell_radius-100.0) {// this is after the insertion, domain can go to negative -100
//                    bound_sum_forces += (cell_radius- x[0]) / x[0] * vdouble2(x[0], 0);
//                }
//                if (x[1] < cell_radius) {
//                    bound_sum_forces += (cell_radius - x[1]) / x[1] * vdouble2(0, x[1]);
//
//                }
//                if (x[1] > length_y - cell_radius - 1) {
//                    bound_sum_forces += (length_y - cell_radius - 1 - x[1]) / x[1] * vdouble2(0, x[1]);
//
//                }
//                if (x[0] > Gamma(length_x - 1) - cell_radius) {
//                    bound_sum_forces += (Gamma(length_x - 1) - cell_radius - x[0]) / x[0] * vdouble2(x[0], 0);
//
//                }
//
//                x = x + (bound_sum_forces); // boundary force
//            }




//            bool free_position = true; // check if the neighbouring position is free
//
//            // if this loop is entered, it means that there is another cell where I want to move
//            for (auto k = euclidean_search(particles.get_query(), x, diameter); k != false; ++k) {
//
//                if (get<id>(*k) !=
//                    get<id>(particles[j])) { // check if it is not the same particle
//                    free_position = false;
//                    break;
//                }
//            }
//
//            // update the position if the place they want to move to is free and not out of bounds
//
//
////            if (freeposition && x[0] > cell_radius && x[0] < Gamma(length_x - 1) && (x[1]) > cell_radius &&
////                (x[1]) < length_y - 1 - cell_radius){
//
//            if (free_position == true) {
//                counttrue = counttrue + 1;
//            }
//            if (free_position == false) {
//                countfalse = countfalse + 1;
//            }






            //get<velocity>(particles[j]) = (sum_forces + random_vector);


            positions(j, 0) = x[0];
            positions(j, 1) = x[1];



        }


        for (int ipos = 0; ipos < particles.size(); ipos++) {

            get<position>(particles)[ipos] = vdouble2(positions(ipos, 0), positions(ipos, 1));

        }




        /* CHECK IF A FOLLOWER DOES NOT BECOME A LEADER
                      * Alternative phenotypic switching if a follower overtakes a leader it becomes a leader and that leader follower.
                      * I will have to be careful when there will be channels because I will have to choose the closest leader
                      * */



        for (int j = 0; j  < particles.size(); j++){
            // find the closest leader


            // so that I would not go through all the cells I will choose the ones that are closer to the front

            // minimum position in x of the leaders

            int min_index = 0;

            for (int i = 1; i < N; ++i) {
                if (get<position>(particles[i])[0] < get<position>(particles[min_index])[0]) {
                    min_index = i;
                }

            }

            // if a follower is eps further in front than the leader, swap their types
            if (get<position>(particles[j])[0] > get<position>(particles[min_index])[0] + eps) {
                // find distance to all the leaders
                double distances[N];
                vdouble2 dist_vector;
                //check which one is the closest
                for (int i = 0; i < N; ++i) {
                    dist_vector = get<position>(particles[j]) - get<position>(particles[i]);
                    distances[i] = dist_vector.norm();

                    int winning_index = 0;
                    for (int i = 1; i < N; ++i) {
                        if (distances[i] < distances[winning_index]) {
                            winning_index = i;
                        }
                    }

                    // if this closest leader is behind that follower, swap them
                    if (get<position>(particles[j])[0] > get<position>(particles[winning_index])[0] + eps) {
                        particle_type::value_type tmp = particles[winning_index];


                        // their position swap

                        vdouble2 temp = get<position>(particles[winning_index]);
                        get<position>(particles[winning_index]) = get<position>(particles[j]);
                        get<position>(particles[j]) = temp;


                    }

                }


            }

        }


// end of phenotype switching














        particles.update_positions(); // not sure if needed here

        if (counter % int(60.0 / dt) == 0) { //  60/dt

        //if (t >1000){
        //if (furthestCell > 980) { // for the one to see how long it takes till they reach the end

        // if (t<5){
        //            //     save at every time step
        #ifdef HAVE_VTK

            // vtkWriteGrid("CellsBehindRepAheadAttreps1p5D5dt0p005MidRemoval", t, particles.get_grid(true)); //.//100sim DATA Parameter Sensitivity Chick Images//
             //          vtkWriteGrid("CellsRepOnlyeps19D12dt0p005MidRemoval", t, particles.get_grid(true)); //.//100sim DATA Parameter Sensitivity Chick Images//
                 vtkWriteGrid("CellsLeaderseps256Foll0D5dt0p005Ablations", t, particles.get_grid(true)); //.//100sim DATA Parameter Sensitivity Chick Images//
                //     vtkWriteGrid("DT0p005Leadersbiasedanddisappearcheck", t, particles.get_grid(true)); //.//100sim DATA Parameter Sensitivity Chick Images//
        #endif
        //}


        // 94, D=1
        // D=6, eps = 19


        }
//////



        // for Mayor's, delete particles greater than 850
        countcellsinarches = 0; // this has to be uncommented!!!
        for (auto p : particles) {
            if (get<position>(p)[0] > Gamma(length_x - 1) - 100.0) {
                countcellsinarches = countcellsinarches + 1;
                // get<alive>(p) = false;
            }
        }


        //} // end of: if (counter % int(60.0 / dt) == 0) {
    }

//    /*
// * return the density of cells in domain_partition parts of the domain
// */
    const int domain_partition = int(round(Gamma(length_x -1 ) / double(100))); // number of intervals of 50 \mu m



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


    return proportions;

}





// parameter analysis
int main() {

    const int number_parameters = 1; // parameter range
    const int sim_num = 1; // it used to be 20 here
    double eps = 256.0; // 0.4, 19.0, 38.0, 56.0, 75.0, 94.0
    //cout << "ignore above" << endl;
//
    //   int num_parts = vector_check_length.size(); // number of parts that I partition my domain
    //  cout << "length " << vector_check_length.size() << endl;
    int num_parts = 11; // for 1080 timesteps
    MatrixXf sum_of_all = MatrixXf::Zero(num_parts, number_parameters); // sum of the values over all simulations



//looping through D
    double D = 5.0;
//        for (int i=1; i < 7; i++) {
//            if (i == 1) {
//                D = 1.0;
//            } else {
//                D = double((i - 1) * 3);
//            }
//            //    cout << "D = " << D << endl;
//
//// looping through beta
    double beta = 0.0;//0.02;
////    for (int i=1; i < 6; i++){
////        beta = 0.01 * double(i);
//            //cout << "beta = " << beta << endl;

    VectorXi vector_check_length =proportions( 0, D, eps, 0.0);; //just to know what the length is

    //    n would correspond to different seeds
    ////     parallel programming
    //WCHAR* folder = L".\\TestFolder";
    double timebutnothing;

    //comment from here

//#pragma omp parallel for
//            for (int n = 0; n < sim_num; n++) {
//
//
//                //initialise the matrix to store the values
//                MatrixXi numbers = MatrixXi::Zero(num_parts, number_parameters);
//
//        //cout << " n = " << n << endl;
//        numbers.block(0, 0, num_parts, 1) = proportions( n, D, eps, beta); // when check proportions
//        //timebutnothing = proportions( n, D, eps, beta); // when check proportions
//
//        // This is what I am using for MATLAB
////        ofstream output2(".//Attr Rep LEAD ONLY Chick Data Files//AttrRepLEADONLYVARYDsepdataChiceps" +to_string(int(eps)) + "D" + to_string(int(D)) + "nvalue" + to_string(n) + ".csv");
////        //
//// DO NOT FORGET TO UNCOMMENT THIS THEN OUTPUT
////        for (int i = 0; i < numbers.rows(); i++) {
////
////            for (int j = 0; j < numbers.cols(); j++) {
////
////              output2 << numbers(i, j) << ", ";
////
////                sum_of_all(i, j) += numbers(i, j);
////
////            }
////            output2 << "\n" << endl;
////        }
//
//    } // end of n simulations
//

    //comment to here




    /*
//    * will store everything in one matrix, the entries will be summed over all simulations
//    */


//    ofstream output3("ChickCiLOnlyeps" +to_string(int(eps)) + "D"+ to_string(int(D)) + ".csv"); // at the end GrowingDomain
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

//}// looping thorugh D or beta
}