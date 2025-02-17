/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

\*---------------------------------------------------------------------------*/

#include "DAInputVolCoord.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAInputVolCoord, 0);
addToRunTimeSelectionTable(DAInput, DAInputVolCoord, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAInputVolCoord::DAInputVolCoord(
    const word inputName,
    const word inputType,
    fvMesh& mesh,
    const DAOption& daOption,
    const DAModel& daModel,
    const DAIndex& daIndex)
    : DAInput(
        inputName,
        inputType,
        mesh,
        daOption,
        daModel,
        daIndex)
{
}

void DAInputVolCoord::run(const scalarList& input)
{
    /*
    Description:
        Assign the input array to the OF's volume coordinates and call movePoints
    */

#ifndef CODI_ADR
    Info << "DAInputVolCoord. " << endl;
    Info << "Setting volume coordinates. " << endl;
#endif

    pointField meshPoints = mesh_.points();

    label counterI = 0;
    forAll(meshPoints, i)
    {
        for (label j = 0; j < 3; j++)
        {
            meshPoints[i][j] = input[counterI];
            counterI++;
        }
    }
    mesh_.movePoints(meshPoints);
    // set this to false to avoid V0 not found error
    // TODO. We need to change this for dynamic mesh cases
    mesh_.moving(false);

    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
