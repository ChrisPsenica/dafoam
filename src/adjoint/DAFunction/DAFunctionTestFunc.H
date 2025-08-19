/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

    Description:
        Child class for mean variable on patches

\*---------------------------------------------------------------------------*/

#ifndef DAFunctionTestFunc_H
#define DAFunctionTestFunc_H

#include "DAFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class DAFunctionTestFunc Declaration
\*---------------------------------------------------------------------------*/

class DAFunctionTestFunc
    : public DAFunction
{

protected:
    /// the area of all total pressure patches
    scalar areaSum_ = -9999.0;

    /// name of the variable
    word varName_;

    /// type of the variable either vector or scalar
    word varType_;

    /// if vector which element/index?
    label index_;

public:
    TypeName("testFunc");
    // Constructors

    //- Construct from components
    DAFunctionTestFunc(
        const fvMesh& mesh,
        const DAOption& daOption,
        const DAModel& daModel,
        const DAIndex& daIndex,
        const word functionName);

    //- Destructor
    virtual ~DAFunctionTestFunc()
    {
    }

    /// calculate the value of objective function
    virtual scalar calcFunction();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
