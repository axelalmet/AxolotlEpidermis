#ifndef TESTTOPDOWNEPIDERMIS_HPP_
#define TESTTOPDOWNEPIDERMIS_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"
#include "CellBasedSimulationArchiver.hpp"


#include "CheckpointArchiveTypes.hpp" // Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "HoneycombMeshGenerator.hpp" //Generates mesh
#include "MutableMesh.hpp"
// #include "FixedRegionPlaneBoundaryCondition.hpp" // Fixed-position boundary condition
#include "PlaneBoundaryCondition.hpp" // Reflective boundary condition
#include "GeneralisedLinearSpringForce.hpp" // General linear spring force that we'll use for steady state
#include "DifferentialWoundAdhesionGeneralisedLinearSpringForce.hpp" // Linear spring force where attraction is weaker than repulsion
#include "PolarityBasedMigrationForce.hpp" // Migration force based on the cell polarity (influenced by mechanics)
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "VoronoiDataWriter.hpp"
#include "VolumeTrackingModifier.hpp"
#include "NodeVelocityWriter.hpp" // Population for cell velocity writer
#include "WildTypeCellMutationState.hpp" // # Default mutation state
#include "TransitCellProliferativeType.hpp" // Default cell proliferation type
#include "ContactInhibitionCellCycleModel.hpp" // Contact-inhibition-based cell cycle
#include "NoCellCycleModel.hpp" // No cell cycle (for when we wound the model)
#include "FakePetscSetup.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "AxolotlEpidermis/TopDown/";
static const double M_DT = 0.005;
static const double M_SS_END_TIME = 50.0;
static const double M_SS_SAMPLING_TIMESTEP = M_SS_END_TIME / M_DT;
static const double M_WOUND_END_TIME = M_SS_END_TIME + 12.0;
static const double M_WOUND_SAMPLING_TIMESTEP = 1.0 / M_DT;

class TestTopDownEpidermis : public AbstractCellBasedTestSuite
{

public:
    void TestWounding()
    {
        /* 
         * First run the simulation to a steady state
         */

        // Define geometry of the model in terms of cells
        unsigned cells_across = 40;
        unsigned cells_up = (unsigned) 160.0 / sqrt(3.0); // This is so that the tissue is actually 800um long

        // Cell cycle parameters for contact inhibition
        double quiescent_fraction = 0.8;
        double equilibrium_area = 3.0*sqrt(3.0)/8.0; // Calculated as area of hexagon with side length 0.5

        // Force parameters
        double spring_stiffness = 30.0;

        HoneycombMeshGenerator generator(cells_across, cells_up, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		//Create shared pointers for cell and mutation states
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

		//Create tissue of cells. Initially we set them all to be differentiated
		std::vector<CellPtr> cells; //Create vector of cells

		for (unsigned i = 0; i < location_indices.size(); i++)
		{
            //Set stochastic duration based cell cycle
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel(); //Don't give them any cell cycle model yet.
            p_cycle_model->SetQuiescentVolumeFraction(quiescent_fraction);
            p_cycle_model->SetEquilibriumVolume(equilibrium_area);
            p_cycle_model->SetDimension(2);
            double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,13) and set in the past
            p_cycle_model->SetBirthTime(-birth_time);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type); //Make cell differentiated
            p_cell->InitialiseCellCycleModel(); // For paranoia really.

            // Set the polarity to be randomly between 0 and 2pi (makes it easier for wounding)
            double polarity = -M_PI + 2.0 * M_PI * RandomNumberGenerator::Instance()->ranf();
            p_cell->GetCellData()->SetItem("polarity", polarity);

            cells.push_back(p_cell);
        }

        // Define cell population
		MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

		//Output data to vtk format so we can visualise it in Paraview
		cell_population.SetWriteVtkAsPoints(true);
		cell_population.AddPopulationWriter<VoronoiDataWriter>();
        cell_population.AddPopulationWriter<NodeVelocityWriter>();

        // Define off-lattice simulation
        OffLatticeSimulation<2> simulator(cell_population);

		//Set output directory
        std::string ss_output_directory = M_OUTPUT_DIRECTORY + "SteadyState";
		simulator.SetOutputDirectory(ss_output_directory);

        // Set timestepping parameters
		simulator.SetDt(M_DT);
		simulator.SetSamplingTimestepMultiple(M_SS_SAMPLING_TIMESTEP); //Sample the simulation at every hour
		simulator.SetEndTime(M_SS_END_TIME); //Hopefully this is long enough for a steady state

        // Add volume-tracking modifier for cell cycle model
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

		//Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
		MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(spring_stiffness);
		p_spring_force->SetCutOffLength(1.5);
		simulator.AddForce(p_spring_force);

		// Define plane boundary conditions for the left, right, top, and bottom row
		c_vector<double, 2> point, normal;

        point(0) = 0.0;
        point(1) = 0.0;
        normal(0) = -1.0;
        normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_left, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_left);

		point(0) = (double)cells_across;
		point(1) = 0.0;
		normal(0) = 1.0;
		normal(1) = 0.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_right, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc_right);

        point(0) = 0.0;
        point(1) = 0.0;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_bottom, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_bottom);

		point(0) = 0.0;
		point(1) = 0.5*sqrt(3.0)*(double)cells_up;
		normal(0) = 0.0;
		normal(1) = 1.0;
		MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_top, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc_top);

		simulator.Solve();

        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator); // Save the simulator
        
        /* Now we wound the model, modelling a loss of cell-cell attraction and an active migration force
         * based on cell polarity, which is influenced by mechanical forces via velocity
         */

        double polarity_remodelling_strength = 5.0; // Allow for sufficiently rapid remodelling to mechanical influences
        std::vector<double> drag_constants = {0.25, 0.5, 0.75, 1.0}
        std::vector<double> migration_force_strengths = {0.0, 1.0, 5.0, 10.0, 50.0};
        std::vector<double> spring_force_multipliers = {0.125, 0.25, 0.5, 1.0};

        for (unsigned i = 0; i < spring_force_multipliers.size(); i++)
        {

            for (unsigned j = 0; j < migration_force_strengths.size(); j++)
            {

                for (unsigned k = 0; k < drag_constants.size(); k++)
                {

                    double adhesion_constant_multiplier = spring_force_multipliers[i];
                    double migration_force_strength = migration_force_strengths[j];
                    double drag_constant = drag_constants[k];

                    OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(ss_output_directory, M_SS_END_TIME);

                    for (AbstractCellPopulation<2>::Iterator cell_iter = p_simulator->rGetCellPopulation().Begin();
                        cell_iter != p_simulator->rGetCellPopulation().End();
                            ++cell_iter)
                    {
                        // Remove the cell cycle by replacing it with a NoCellCycleModel
                                    //Set stochastic duration based cell cycle
                        NoCellCycleModel * p_cycle_model = new NoCellCycleModel(); //Don't give them any cell cycle model yet.
                        p_cycle_model->SetDimension(2);
                        cell_iter->SetCellCycleModel(p_cycle_model);
                        cell_iter->InitialiseCellCycleModel(); // For paranoia really.

                        unsigned node_index = p_simulator->rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter);

                        double y = p_simulator->rGetCellPopulation().rGetMesh().GetNode(node_index)->rGetLocation()[1];
                        
                        // We need to initialise the applied for contribution for the writer, I think
                        Node<2>* p_node = p_simulator->rGetCellPopulation().rGetMesh().GetNode(node_index);

                        p_node->AddAppliedForceContribution(zero_vector<double>(2));

                        // Kill the cell if it's above the upper half of the boundary
                        if (y > 0.25*sqrt(3.0)*(double)cells_up)
                        {
                            cell_iter->Kill();
                        }
                    }

                    p_simulator->rGetCellPopulation().RemoveDeadCells();
                    p_simulator->rGetCellPopulation().Update();

                    // Change the drag constant now
                    p_simulator->rGetCellPopulation().SetDampingConstantNormal(drag_constant);

                    // Add the new force
                    p_simulator->RemoveAllForces();

                    // Add spring force to account for loss of adhesion in wounding
                    MAKE_PTR(DifferentialWoundAdhesionGeneralisedLinearSpringForce<2>, p_wound_spring_force);
                    p_wound_spring_force->SetMeinekeSpringStiffness(spring_stiffness);
                    p_wound_spring_force->SetAdhesionSpringConstantMultiplier(adhesion_constant_multiplier);
                    p_wound_spring_force->SetCutOffLength(1.5);
                    p_simulator->AddForce(p_wound_spring_force);

                    // Add polarity-based migration force
                    MAKE_PTR(PolarityBasedMigrationForce<2>, p_migration_force);
                    p_migration_force->SetMigrationForceStrength(migration_force_strength);
                    p_migration_force->SetPolarityRemodellingStrength(polarity_remodelling_strength);
                    simulator.AddForce(p_migration_force);

                    //Set output directory
                    std::stringstream out;
                    out << "Wounding/ADHESION_" << adhesion_constant_multiplier << "_MIGRATION_" << migration_force_strength << "_DRAG_" << drag_constant;
                    std::string wound_output_directory = M_OUTPUT_DIRECTORY + out.str();

                    p_simulator->SetOutputDirectory(wound_output_directory);

                    p_simulator->SetEndTime(M_WOUND_END_TIME);
                    p_simulator->SetSamplingTimestepMultiple(M_WOUND_SAMPLING_TIMESTEP);

                    // Run the wounding part
                    p_simulator->Solve();

                    delete p_simulator;

                    //Tidying up
                    SimulationTime::Instance()->Destroy();
                    SimulationTime::Instance()->SetStartTime(M_SS_END_TIME);
                }
            }
        }
    }

};

#endif /* TESTTOPDOWNEPIDERMIS_HPP_ */
