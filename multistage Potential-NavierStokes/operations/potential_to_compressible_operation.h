//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//
//

#ifndef KRATOS_POTENTIAL_TO_COMPRESSIBLE_OPERATION_INCLUDED
#define KRATOS_POTENTIAL_TO_COMPRESSIBLE_OPERATION_INCLUDED


// System includes


// External includes


// Project includes
#include "operations/operation.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @author Miguel Maso Sotomayor
 * @class PotentialToCompressibleOperation
 * @brief Copy the properties from one model part to another.
 * @details The properties of the elements and conditions of the destination model part are replaced by the copies.
 */
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) PotentialToCompressibleOperation : public Operation
{
public:
    ///@name Pointer Definition
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PotentialToCompressibleOperation);

    ///@}
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    ///@brief Default constructor.
    PotentialToCompressibleOperation() : Operation() {}

    /**
     * @brief Constructor with Model and Parameters
     * @param rModel Reference of the Model
     * @param ModelerParameters Parameters of the Modeler
     */
    PotentialToCompressibleOperation(
        Model& rModel,
        Parameters ModelerParameters);

    /**
     * @brief Constructor with origin and destination model parts
     * @param rOriginModelPart Reference of the origin ModelPart
     * @param rDestinationModelPart Reference of the destination ModelPart
     */
    PotentialToCompressibleOperation(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart);

    ///@brief Copy constructor.
    PotentialToCompressibleOperation(PotentialToCompressibleOperation const& rOther) = delete;

    ///@brief Destructor.
    ~PotentialToCompressibleOperation() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@brief Assignment operator.
    PotentialToCompressibleOperation & operator=(PotentialToCompressibleOperation const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates the Modeler Pointer and returns a pointer to a new PotentialToCompressibleOperation
     * @param rModel Reference of the Model
     * @param ModelerParameters Parameters of the discretization
     * @return a Pointer to the new Modeler
     */
    Modeler::Pointer Create(
        Model& rModel,
        const Parameters ModelParameters) const override;

    /**
     * @brief Get the Default Parameters object
     * @return * This const 
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Stages
    ///@{

    /**
     * @brief Create a copy of the properties from the origin ModelPart to the destination ModelPart
     */
    void Excecute() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "PotentialToCompressibleOperation";
    }

    ///@}

private:

    ///@name Member Variables
    ///@{

    Model* mpModel = nullptr;

    ///@}
    ///@name Operations
    ///@{

    void RecursivelyCopyProperties(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart);

    template<class TContainerType>
    void ReplaceProperties(
        TContainerType& rContainer,
        const ModelPart& rModelPart);

    ///@}
};

///@}
///@name Input and output
///@{

///@brief input stream function
inline std::istream& operator>>(std::istream& rIStream,
    PotentialToCompressibleOperation& rThis)
{
    return rIStream;
}

///@brief output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const PotentialToCompressibleOperation& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

} // namespace Kratos.

#endif //KRATOS_POTENTIAL_TO_COMPRESSIBLE_OPERATION_INCLUDED defined
