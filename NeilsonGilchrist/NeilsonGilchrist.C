/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "NeilsonGilchrist.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::NeilsonGilchrist<CloudType>::applyToPatch
(
    const label globalPatchi
) const
{
    forAll(patchIDs_, i)
    {
        if (patchIDs_[i] == globalPatchi)
        {
            return i;
        }
    }

    return -1;
}


template<class CloudType>
void Foam::NeilsonGilchrist<CloudType>::write()
{
    if (NeilsonGilchristPtr_.valid())
    {
        NeilsonGilchristPtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "NeilsonGilchristPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NeilsonGilchrist<CloudType>::NeilsonGilchrist
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    NeilsonGilchristPtr_(nullptr),
    patchIDs_(),
    n_(this->coeffDict().template lookupOrDefault<scalar>("n", 5)),
    epsilonN_(this->coeffDict().template lookupOrDefault<scalar>("epsilonN", 7.85e10)),
    rhoN_(this->coeffDict().template lookupOrDefault<scalar>("rhoN", 1.77e10)),
    alpha0_(this->coeffDict().template lookupOrDefault<scalar>("alpha0", 0.3)),
    K_(this->coeffDict().template lookupOrDefault<scalar>("K", 0))

{
    const wordList allPatchNames = owner.mesh().boundaryMesh().names();
    wordList patchName(this->coeffDict().lookup("patches"));

    labelHashSet uniquePatchIDs;
    forAllReverse(patchName, i)
    {
        labelList patchIDs = findStrings(patchName[i], allPatchNames);

        if (patchIDs.empty())
        {
            WarningInFunction
                << "Cannot find any patch names matching " << patchName[i]
                << endl;
        }

        uniquePatchIDs.insert(patchIDs);
    }

    patchIDs_ = uniquePatchIDs.toc();

    // Trigger creation of the NG field
    preEvolve();
}


template<class CloudType>
Foam::NeilsonGilchrist<CloudType>::NeilsonGilchrist
(
    const NeilsonGilchrist<CloudType>& pe
)
:
    CloudFunctionObject<CloudType>(pe),
    NeilsonGilchristPtr_(nullptr),
    patchIDs_(pe.patchIDs_),
    n_(pe.n_),
    epsilonN_(pe.epsilonN_),
    rhoN_(pe.rhoN_),
    alpha0_(pe.alpha0_),
    K_(pe.K_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NeilsonGilchrist<CloudType>::~NeilsonGilchrist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::NeilsonGilchrist<CloudType>::preEvolve()
{
    if (NeilsonGilchristPtr_.valid())
    {
        NeilsonGilchristPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        NeilsonGilchristPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "NG",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimVolume, 0.0)
            )
        );
    }
}


template<class CloudType>
void Foam::NeilsonGilchrist<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp,
    bool&
)
{
    const label patchi = pp.index();

    const label localPatchi = applyToPatch(patchi);

    if (localPatchi != -1)
    {
        vector nw;
        vector Up;

        // patch-normal direction
        this->owner().patchData(p, pp, nw, Up);

        // particle velocity relative to patch
        const vector& U = p.U() - Up;

        // quick reject if particle travelling away from the patch
        if ((nw & U) < 0)
        {
            return;
        }

        const scalar magU = mag(U);
        const vector Udir = U/magU;

        // determine impact angle, and constants
        const scalar alpha = mathematical::pi/2.0 - acos(nw & Udir);

 	const scalar M = p.nParticle()*p.mass();

        const label patchFacei = pp.whichFace(p.face());
        scalar& NG = NeilsonGilchristPtr_->boundaryFieldRef()[patchi][patchFacei];

	if (alpha < alpha0_)
        {
		NG += ((0.5 * M * pow(magU,2) * pow(cos(alpha),2) * sin(n_*alpha))/(rhoN_)) + ((0.5 * M * pow((magU * sin(alpha) - K_),2))/(epsilonN_));
	}
        else 
        { 
		NG += ((0.5 * M * pow(magU,2) * pow(cos(alpha),2))/(rhoN_)) + ((0.5 * M * pow((magU * sin(alpha) - K_),2))/(epsilonN_)); 
	}
    }
}


// ************************************************************************* //
