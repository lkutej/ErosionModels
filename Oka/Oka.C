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

#include "Oka.H"

// * * * * * * * * * * * * * Protectd Member Functions * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::Oka<CloudType>::applyToPatch
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
void Foam::Oka<CloudType>::write()
{
    if (OkaPtr_.valid())
    {
        OkaPtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "OkaPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Oka<CloudType>::Oka
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    OkaPtr_(nullptr),
    patchIDs_(),
    n1_(readScalar(this->coeffDict().lookup("n1"))),
    n2_(readScalar(this->coeffDict().lookup("n2"))),
    Hv_(readScalar(this->coeffDict().lookup("Hv"))),
    K_(readScalar(this->coeffDict().lookup("K"))),
    k1_(readScalar(this->coeffDict().lookup("k1"))),
    k2_(readScalar(this->coeffDict().lookup("k2"))),
    k3_(readScalar(this->coeffDict().lookup("k3")))
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
Foam::Oka<CloudType>::Oka
(
    const Oka<CloudType>& pe
)
:
    CloudFunctionObject<CloudType>(pe),
    OkaPtr_(nullptr),
    patchIDs_(pe.patchIDs_),
    n1_(pe.n1_),
    n2_(pe.n2_),
    Hv_(pe.Hv_),
    K_(pe.K_),
    k1_(pe.k1_),
    k2_(pe.k2_),
    k3_(pe.k3_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::Oka<CloudType>::~Oka()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::Oka<CloudType>::preEvolve()
{
    if (OkaPtr_.valid())
    {
        OkaPtr_->primitiveFieldRef() = 0.0;
    }
    else
    {
        const fvMesh& mesh = this->owner().mesh();

        OkaPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    this->owner().name() + "O",
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
void Foam::Oka<CloudType>::postPatch
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
        scalar& O = OkaPtr_->boundaryFieldRef()[patchi][patchFacei];

        const scalar E90 = K_ * pow(Hv_,k1_)*pow(magU,k2_)*pow(p.d(),k3_);

        const scalar FuncG = pow(sin(alpha),n1_) * pow(1 + Hv_ *(1- sin(alpha)),n2_);
	    
        O += M * FuncG * E90;

    }
}


// ************************************************************************* //
