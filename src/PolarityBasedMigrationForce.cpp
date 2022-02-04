/*

Copyright (c) 2005-2020, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "PolarityBasedMigrationForce.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
PolarityBasedMigrationForce<DIM>::PolarityBasedMigrationForce()
    : AbstractForce<DIM>(),
    mMigrationForceStrength(DOUBLE_UNSET),
    mPolarityRemodellingStrength(DOUBLE_UNSET)
{
    if (DIM != 2)
    {
        EXCEPTION("This force is only implemented in 2D.");
    }
}

template<unsigned DIM>
PolarityBasedMigrationForce<DIM>::~PolarityBasedMigrationForce()
{
}

// Method to get the strength of the migration force
template<unsigned DIM>
double PolarityBasedMigrationForce<DIM>::GetMigrationForceStrength()
{
    return mMigrationForceStrength;
}

// Method to set the strength of migration force
template<unsigned DIM>
void PolarityBasedMigrationForce<DIM>::SetMigrationForceStrength(double migrationForceStrength)
{
    mMigrationForceStrength = migrationForceStrength;
}

// Method to get the strength of the migration force
template<unsigned DIM>
double PolarityBasedMigrationForce<DIM>::GetPolarityRemodellingStrength()
{
    return mPolarityRemodellingStrength;
}

// Method to set the strength of migration force
template<unsigned DIM>
void PolarityBasedMigrationForce<DIM>::SetPolarityRemodellingStrength(double polarityRemodellingStrength)
{
    mPolarityRemodellingStrength = polarityRemodellingStrength;
}


template<unsigned DIM>
double PolarityBasedMigrationForce<DIM>::UpdateCellPolarity(AbstractCellPopulation<DIM>& rCellPopulation, unsigned nodeIndex, double cellPolarity)
{
    double new_cell_polarity;

    // This really only makes sense if we're working with a centre-based population
    if (bool(dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation)))
    {
        double damping_constant = static_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation)->GetDampingConstant(nodeIndex);

        // Get the time step
        double dt = SimulationTime::Instance()->GetTimeStep(); 

        // Get the applied force to work out the velocity
        c_vector<double, DIM> velocity = rCellPopulation.rGetMesh().GetNode(nodeIndex)->rGetAppliedForce() / damping_constant;

        // Get the angle the velocity makes with respect to the x-axis
        double velocity_angle = atan2(velocity[1], velocity[0]); // Use atan2 to adjust the angle between -pi and pi

        // Return the new angle
        new_cell_polarity =  cellPolarity + dt * mPolarityRemodellingStrength * (velocity_angle - cellPolarity);
        
    }
    else // If not, we just don't update it
    {
        new_cell_polarity = cellPolarity;
    }

    return new_cell_polarity;

}

template<unsigned DIM>
void PolarityBasedMigrationForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

        // Make sure we don't do this to dead cells
        if (!cell_iter->IsDead())
        {
            // Get the node index
            unsigned current_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            // Get the pointer to the node
            Node<DIM>* p_node = rCellPopulation.GetNode(current_index);

            double old_polarity = cell_iter->GetCellData()->GetItem("polarity");

            double new_polarity = UpdateCellPolarity(rCellPopulation, current_index, old_polarity);

            // Set the new cell polarity
            cell_iter->GetCellData()->SetItem("polarity", new_polarity);

            // Add the migration force (only implemented in 2D, lol)
            c_vector<double, 2> migration_force = zero_vector<double>(2);
            
            migration_force[0] = mMigrationForceStrength * cos(new_polarity);
            migration_force[1] = mMigrationForceStrength * sin(new_polarity);

            p_node->AddAppliedForceContribution(migration_force);
        }
            
    }
}

template<unsigned DIM>
void PolarityBasedMigrationForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include
	*rParamsFile <<  "\t\t\t<MigrationForceStrength>"<<  mMigrationForceStrength << "</MigrationForceStrength> \n" ;
	*rParamsFile <<  "\t\t\t<PolarityRemodellingStrength>"<<  mPolarityRemodellingStrength << "</PolarityRemodellingStrength> \n" ;

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class PolarityBasedMigrationForce<1>;
template class PolarityBasedMigrationForce<2>;
template class PolarityBasedMigrationForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PolarityBasedMigrationForce)
