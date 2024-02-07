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

#include "FinnieExtended.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::FinnieExtended<CloudType>::applyToPatch
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
void Foam::FinnieExtended<CloudType>::write()
{
    if (FinnieExtendedPtr_.valid())
    {
        FinnieExtendedPtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "FinnieExtendedPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FinnieExtended<CloudType>::FinnieExtended
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    FinnieExtendedPtr_(nullptr),
    patchIDs_(),
    P_(readScalar(this->coeffDict().lookup("P"))),
    K_(this->coeffDict().template lookupOrDefault<scalar>("K", 1)), //sphere: 1 square: 0.4
    epsilonA_(this->coeffDict().template lookupOrDefault<scalar>("epsilonA", 7.85e10)), 
    Fs_(this->coeffDict().template lookupOrDefault<scalar>("Fs", 1)), //fully round: 0.25 fully sharp: 1    
    C_(this->coeffDict().template lookupOrDefault<scalar>("C", 0.0042)),
    Utsh_(this->coeffDict().template lookupOrDefault<scalar>("Utsh", 7.4))	

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

    // Trigger creation of the A field
    preEvolve();
}


template<class CloudType>
Foam::FinnieExtended<CloudType>::FinnieExtended
(
    const FinnieExtended<CloudType>& pe
)
:
    CloudFunctionObject<CloudType>(pe),
    FinnieExtendedPtr_(nullptr),
    patchIDs_(pe.patchIDs_),
    P_(pe.P_),
    K_(pe.K_),
    epsilonA_(pe.epsilonA_),
    Fs_(pe.Fs_),
    C_(pe.C_),
    Utsh_(pe.Utsh_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FinnieExtended<CloudType>::~FinnieExtended()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::FinnieExtended<CloudType>::preEvolve()
{
    if (FinnieExtendedPtr_.valid())
    {
        FinnieExtendedPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        FinnieExtendedPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "FE",
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
void Foam::FinnieExtended<CloudType>::postPatch
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
        scalar& FE = FinnieExtendedPtr_->boundaryFieldRef()[patchi][patchFacei];

	if (alpha < atan(K_) && magU*sin(alpha) > Utsh_)
        {
		FE += Fs_ * (((0.5*M*pow((magU*sin(alpha) - Utsh_),2))/(epsilonA_)) + (C_*((M*pow(magU,2)*sin(alpha)*(2*K_*cos(alpha)-sin(alpha)))/(2*pow(K_,2)*P_))));
	}
        else if (alpha >= atan(K_) && magU*sin(alpha) > Utsh_)
        { 
		FE += Fs_ * (((0.5*M*pow((magU*sin(alpha) - Utsh_),2))/(epsilonA_)) + (C_*((M * pow(magU,2)*pow(cos(alpha),2))/(2*P_)))); 
	}
	else
	{
		FE += 0;
	}

    }
}


// ************************************************************************* //
