#include "sequenceFactory.h"

#include "elements/updateForcesElement.h"
#include "elements/ouoperatorElement.h"
#include "elements/constrainElement.h"
#include "elements/updatePositionsElement.h"
#include "elements/updateVelocitiesElement.h"
#include "elements/computeGlobalsElement.h"
#include "elements/writeStateElement.h"
#include "elements/writeEnergyElement.h"
#include "elements/helloWorldElement.h"
#include "elements/saveSimulationStateElement.h"
#include "elements/saveEnergyElement.h"
#include "elements/restoreSimulationStateElement.h"
#include "elements/metropolisAcceptanceTestElement.h"
#include "elements/drawVelocitiesFromGaussianElement.h"

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/inputrec.h"

/*
    method that creates the PreRunSequence,
    ***Currently***: Presequence is created depending on the information
    given by ir (mdp file) and bConstrain

    ***In Future***: A Seperate file where the user lists the elements to construct an integrator
    will be used to create the preSequence.
*/
std::list<Element*> SequenceFactory::createPreRunSequence(t_inputrec *ir, bool bConstrain)
{
    std::list<Element*> sequence;

    if (bConstrain)
    {
        bool bHalfStep=FALSE;
        if (EI_CUSTOM(ir->eI) && ir->eCT == ectLF_NVE)
        {
            bHalfStep = TRUE;
        }
        else if (EI_CUSTOM(ir->eI) && ir->eCT == ectVV_NVE)
        {
            bHalfStep = FALSE;
        }
        else if (EI_CUSTOM(ir->eI) && ir->eCT == ectLANGEVIN)
        {
            bHalfStep = FALSE;
        }
        else if (EI_CUSTOM(ir->eI) && ir->eCT == ectLF_HMC)
        {
            bHalfStep = TRUE;
        }
        sequence.emplace_back(new InitialConstrain(bHalfStep));
    }
    return sequence;
}

/*
    method that creates the runSequence, which is the algorithmic
    representation of the steps that get translated to the code in the
    integrator object.

    ***Currenty***: runsequence is created depending on the information
    given by ir (mdp file) and bConstrain

    in the mdp file the "integrator" field has an option "custom"
    another field named "custom-type" has been added to the mdfile, 
    which determines the custom integrator type.
    Using the information from "custom-type", a Seqence of element is constrcuted.

    ***In Future***: Instead of giving the custom-integrator type, the user will be able to list
    the elements and the arguments in a file. This file will be parsed, and the information from that file
    will be used to construct the sequence of elements.
*/
std::list<Element*> SequenceFactory::createRunSequence(t_inputrec *ir, bool bConstrain)
{
    std::list<Element*> sequence;
    
    if (EI_CUSTOM(ir->eI) && ir->eCT == ectLF_NVE)
    {
        /*
        * A LeapFrog update Algorithm, It also applies contrain, if constrain is set.
        * STATE -> v(t-1/2dt), x(t), f(t)
        * v(t+1/2dt) = v(t-1/2dt) + (1/m)f(t)dt
        * x(t+dt) = x(t) + v(t+1/2dt)dt
        * STATE -> v(t+1/2dt), x(t+dt), f(t+dt)
        */
        
        /* 
            ***!!! Temporary Fix !!!***
            The following booleans are needed to tell computeglobals what
            global variables to calculate.
            The logic is the same as in the standard gromacs compute_globals function

            Future Work: Make the logic simpler for the user. 
        */

        gmx_bool bTemp=TRUE;        // Calculate the kinetic energy and temperature
        gmx_bool bScaleEkin=FALSE;  // scale the kinetic energy with the nose-hoover chain scalar
        gmx_bool bReadEkin=FALSE;   // bReadEkin=TRUE will make bEkinAveVel=TRUE
        
        /* 
            if bEkinAveVel is FALSE:
              the Kinetic energy is calculated via the average of the halfstep kinetic energies
            if bEkinAveVel is TRUE:
              the kinetic energy is calculated from the given velocities
            bEkinAveVel is True if either bReadEkin=True, or the relevant cglo_flag is true
        */
        gmx_bool bEkinAveVel=FALSE; 

        gmx_bool bEnergy=FALSE;     // Calculate dispersion correction and add it to the potential energy
        gmx_bool bStopCM=FALSE;     // Remove com motion from velocities
        gmx_bool bPres=TRUE;        // Calculate dispersion correction and pressure
        gmx_bool bConstraint=FALSE; // Does the same thing as bPres (May not be neccessary)
        bool bInterSimSignal=false; 
        bool bIntraSimSignal=false;
         
        sequence.emplace_back(new UpdateForces());
        sequence.emplace_back(new WriteState());
        sequence.emplace_back(new UpdateVelocities(1.0));
        sequence.emplace_back(new UpdatePositions(1.0));

        if (bConstrain)
        {
            // Constrain the positions
            // it also constraints the half step velocities
            /*
                ***!!! Note!!!***
                Ideally we want to have two seperate elements where one does 
                constraints only the velocities and the other only the positions. 
                Such as:
                * ConstrainVelocity()
                * ConstrainPosition()
                
                Currently we have an element called "ConstrainPositionAndVelocity"
                that does position constraint and also constrain the halfstep velocities
                using the information of x(t) and x(t+dt). 
                Constraining the halfstep velocities are done in a similar way as the 
                position constrain. Thus the algorithms is tied into the position 
                constrain algorithm, which makes it non trivial to seperate it.

                A further discussian is needed about seperating the velocity constrain
                from the position constrain without dublicating code
            */
            sequence.emplace_back(new ConstrainPositionAndVelocity());
        }
        sequence.emplace_back(new ComputeGlobals(bEnergy, bReadEkin, bScaleEkin,
                                                 bEkinAveVel, bStopCM, bTemp, bPres,
                                                 bConstraint, bInterSimSignal, bIntraSimSignal));
        sequence.emplace_back(new WriteEnergy());
    }
    else if (EI_CUSTOM(ir->eI) && ir->eCT == ectVV_NVE)
    {
        /*
        * Velocity Verlet update Algorithm, , It also applies contrain, if constrain is set.
        * -> STATE v(t), x(t), f(t)
        * v(t+1/2dt) = v(t) + (1/m)f(t)dt
        * x(t+dt) = x(t) + v(t+1/2dt)dt
        * v(t+dt) = v(t+1/2dt) + (1/m)f(t+dt)dt
        * -> STATE v(t+dt), x(t+dt), f(t+dt)
        */

        // See explanation above
        gmx_bool bTemp=TRUE;
        gmx_bool bScaleEkin=FALSE;  
        gmx_bool bReadEkin=FALSE;   
        gmx_bool bEkinAveVel=TRUE; 
        gmx_bool bEnergy=FALSE;
        gmx_bool bStopCM=FALSE;
        gmx_bool bPres=TRUE;
        gmx_bool bConstraint=FALSE;
        bool bInterSimSignal=false;
        bool bIntraSimSignal=false;

        // runs only at first step to have forces at t = 0
        sequence.emplace_back(new UpdateForces(0, true)); 
        sequence.emplace_back(new ComputeGlobals(bEnergy, bReadEkin, bScaleEkin,
                                                 bEkinAveVel, bStopCM, bTemp, bPres,
                                                 bConstraint, bInterSimSignal, bIntraSimSignal));
        sequence.emplace_back(new WriteEnergy());
        sequence.emplace_back(new WriteState());
        sequence.emplace_back(new UpdateVelocities(0.5));
        sequence.emplace_back(new UpdatePositions(1.0));
        
        if (bConstrain)
        {
            // Constraint the position
            // it also constraint the half step velocities
            sequence.emplace_back(new ConstrainPositionAndVelocity());
        }
        
        sequence.emplace_back(new UpdateForces());
        sequence.emplace_back(new UpdateVelocities(0.5));
        
        if (bConstrain)
        {
            // This constrains the full step velocity
            // using x(t+dt) and unconstrained v(t+dt)
            sequence.emplace_back(new ConstrainVelocity());
        }
    }
    else if (EI_CUSTOM(ir->eI) && ir->eCT == ectLANGEVIN)
    {
        /*
         * Langevin Integrator algorithm as in:
         * D. A. Sivak et. al. ,J. Phys. Chem. B, 118, 6466-6474, 2014
         * -> STATE v(t), x(t), f(t)
         * v(t+1/4dt) = a*v(t) + sqrt[(1-a^2)/(beta*m)]*rand(t)
         * v(t+1/2dt) = v(t+1/4dt) + (b*dt)/2 * f(t)/m
         * x(t+dt)   = x(t) + b*dt * v(t+1/2dt)
         * f(t+dt)
         * v(t+3/4dt) = v(t+1/2) + (b*dt)/2 * f(t+dt)/m
         * v(t+dt)   = sqrt[a]*v(t+3/4dt) + sqrt[(1-a)/(beta*m)]*rand(t+dt)
         * Cache STATE
         * -> STATE v(t+dt), x(t+dt), f(t+dt)
         *
         * a = exp[-gamma * dt / 2]
         * b = 1, for now
        */

        gmx_bool bTemp=TRUE;
        gmx_bool bScaleEkin=FALSE;  
        gmx_bool bReadEkin=FALSE; 
        gmx_bool bEkinAveVel=FALSE;
        gmx_bool bEnergy=FALSE;
        gmx_bool bStopCM=FALSE;
        gmx_bool bPres=TRUE;
        gmx_bool bConstraint=FALSE;
        bool bInterSimSignal=false;
        bool bIntraSimSignal=false;

        /*
         * In the current Langevin code, these lines were executed before
         * every update. Do we have them somewhere now? TODO: CHECK!
         */
        /* We need to update the NMR restraint history when time averaging is used * /
        if (state->flags & (1<<estDISRE_RM3TAV))
        {
            update_disres_history(fcd, &state->hist);
        }
        if (state->flags & (1<<estORIRE_DTAV))
        {
            update_orires_history(fcd, &state->hist);
        }
        */

        // runs only at first step to have forces at t = 0
        sequence.emplace_back(new UpdateForces(0, true));
        sequence.emplace_back(new ComputeGlobals(bEnergy, bReadEkin, bScaleEkin,
                                                 bEkinAveVel, bStopCM, bTemp, bPres,
                                                 bConstraint, bInterSimSignal, bIntraSimSignal));
        sequence.emplace_back(new WriteEnergy());
        sequence.emplace_back(new WriteState());

        sequence.emplace_back(new OUOperator(0.5));
        sequence.emplace_back(new UpdateVelocities(0.5));
        sequence.emplace_back(new UpdatePositions(1.0));
        sequence.emplace_back(new UpdateForces());
        sequence.emplace_back(new UpdateVelocities(0.5));
        sequence.emplace_back(new OUOperator(0.5));
    }
    else if (EI_CUSTOM(ir->eI) && ir->eCT == ectLF_HMC)
    {
        /*
        * 10 LeapFrog updates in one step, It also applies contrain, if constrain is set.
        * for i in MDStep:
        *   STATE -> v(t-1/2dt), x(t), f(t)
        *   v(t+1/2dt) = v(t-1/2dt) + (1/m)f(t)dt
        *   x(t+dt) = x(t) + v(t+1/2dt)dt
        *   STATE -> v(t+1/2dt), x(t+dt), f(t+dt)
        */
        
        int MDStep = 10; // Number of MDStep per HMC step

        gmx_bool bTemp=TRUE;        
        gmx_bool bScaleEkin=FALSE;  
        gmx_bool bReadEkin=FALSE;   
        gmx_bool bEkinAveVel=FALSE; 
        gmx_bool bEnergy=FALSE;     
        gmx_bool bStopCM=FALSE;    
        gmx_bool bPres=TRUE;       
        gmx_bool bConstraint=FALSE; 
        bool bInterSimSignal=false; 
        bool bIntraSimSignal=false;
        
        // Randomize Velocities
        sequence.emplace_back(new DrawVelocitiesFromGaussian());
        /* 
        The compute global is called to calculate ekinh, using the newly generated velocities.
        after another call to compute global after a velocity update, the ekinh will be copied
        to ekinh_old and ekinh will be recalculated with the new velocities.
        Thus when the temperature of kinetic energies are calculated it will use the average of
        ekinh and ekinh_old

        */
        sequence.emplace_back(new ComputeGlobals(bEnergy, bReadEkin, bScaleEkin,
                                                 bEkinAveVel, bStopCM, bTemp, bPres,
                                                 bConstraint, bInterSimSignal, bIntraSimSignal));
        /* 
           ***!Note:***!
           We might need to constrain the velocities after drawing 
           the velocities from a Gaussian distribution

           However since at this stage the velocities are the halfstep velocities,
           v(t-1/2dt) to constrain at this stage we need,
           x(t-dt) and x(t). and apply the ConstrainPositionAndVelocity() Element Afterwards
           
           In principle we could use the same code construct as in the intialConstrain Element,
           But it might not be ideal to calculate the x(t-dt) at every HMC step.

           Currently the alternative would be to store x(t-dt) and x(t) when saveing and restoring the
           simulation state, which would get rid of explicitly calculating x(t-dt)

           This problem would only occur when using a leapfrog, it should be fine for a velocity verlet 
           integrator.
        */
        sequence.emplace_back(new SaveSimulationState(1,"a"));
        // The loop for the MD iteration. This loop is done every HMC step
        for(int i=0; i<= MDStep; i++)
        {
            if (i == MDStep)
            {
                // Save the state at the last step t= MDStep
                sequence.emplace_back(new SaveSimulationState(1,"b"));
            }
            sequence.emplace_back(new UpdateForces());
            sequence.emplace_back(new WriteState());
            sequence.emplace_back(new UpdateVelocities(1.0));
            sequence.emplace_back(new UpdatePositions(1.0));

            if (bConstrain)
            {
                sequence.emplace_back(new ConstrainPositionAndVelocity());
            }
            sequence.emplace_back(new ComputeGlobals(bEnergy, bReadEkin, bScaleEkin,
                                                     bEkinAveVel, bStopCM, bTemp, bPres,
                                                     bConstraint, bInterSimSignal, bIntraSimSignal));
            if (i == 0)
            {
                // energies of state t=0
                sequence.emplace_back(new SaveEnergy(1,"a"));
            }
            sequence.emplace_back(new WriteEnergy());
        }
        /* Store the current state for the acceptanceTest to use
        Note: If this will be always the current state, than there might not be the need to store it.
        SIn that case the Acceptance test would always assume the state to compare the old state
        is always the current state.
        */
        // Save energies at t=MDStep
        sequence.emplace_back(new SaveEnergy(1,"b"));
        // Restore old state if test fails
        sequence.emplace_back(new MetropolisAcceptanceTest("a", "b"));
    }
    return sequence;
}

/*
    Method that creates the PostRunSequence,
    Was added for potential future purposes, currently returns empty sequence
*/
std::list<Element*> SequenceFactory::createPostRunSequence(t_inputrec *ir,
                                                           bool bConstrain)
{
    std::list<Element*> sequence;

    // add elements here

    return sequence;
}
