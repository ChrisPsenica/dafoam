/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v4

\*---------------------------------------------------------------------------*/

#include "DAInputPatchVx.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAInputPatchVx, 0);
addToRunTimeSelectionTable(DAInput, DAInputPatchVx, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAInputPatchVx::DAInputPatchVx(
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

void DAInputPatchVx::run(const scalarList& input)
{
    /*
    Description:
        Assign the input array to OF's state variables
    */

    // NOTE: we need to first update DAGlobalVar::patchVelocity here, so that daFunction-force
    // can use it to compute the force direction.
    wordList patchNames;
    dictionary patchVSubDict = daOption_.getAllOptions().subDict("inputInfo").subDict(inputName_);
    patchVSubDict.readEntry<wordList>("patches", patchNames);

    scalar Vx = input[0];

    volVectorField& U = const_cast<volVectorField&>(mesh_.thisDb().lookupObject<volVectorField>("U"));

    // now we assign UxNew and UyNew to the U patches
    forAll(patchNames, idxI)
    {
        word patchName = patchNames[idxI];
        label patchI = mesh_.boundaryMesh().findPatchID(patchName);
        if (mesh_.boundaryMesh()[patchI].size() > 0)
        {
            if (U.boundaryField()[patchI].type() == "fixedValue")
            {
                forAll(U.boundaryField()[patchI], faceI)
                {
                    U.boundaryFieldRef()[patchI][faceI][0] = Vx;
                }
            }
            else if (U.boundaryField()[patchI].type() == "inletOutlet")
            {
                mixedFvPatchField<vector>& inletOutletPatch =
                    refCast<mixedFvPatchField<vector>>(U.boundaryFieldRef()[patchI]);

                forAll(U.boundaryField()[patchI], faceI)
                {
                    inletOutletPatch.refValue()[faceI][0] = Vx;
                }
            }
            else
            {
                FatalErrorIn("DAInputPatchVx::run")
                    << "patch type not valid! only support fixedValue or inletOutlet"
                    << exit(FatalError);
            }
        }
    }
    U.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
