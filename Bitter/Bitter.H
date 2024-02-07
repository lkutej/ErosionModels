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

Class
    Foam::Bitter

Description
    Creates particle erosion field, W

SourceFiles
    Bitter.C

\*---------------------------------------------------------------------------*/

#ifndef Bitter_H
#define Bitter_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleErosion Declaration
\*---------------------------------------------------------------------------*/
template<class CloudType>
class Bitter
:
    public CloudFunctionObject<CloudType>
{
    // Private Data

        // Typedefs

            //- Convenience typedef for parcel type
            typedef typename CloudType::parcelType parcelType;


        //- Particle erosion field
        autoPtr<volScalarField> BitterPtr_;

        //- List of patch indices to post-process
        labelList patchIDs_;

        //- elastic load limit
        scalar y_;

        //- deformation wear factor
        scalar epsilonB_;

        //- density
        scalar d_;

        //- cutting wear factor
        scalar rhoB_;

	//- decision angle
        scalar alpha0_;

	//- poission ratio particles
        scalar q1_;

	//- poission ratio body
        scalar q2_;

	//- Young modulus particles
        scalar E1_;

	//- Young modulus body
        scalar E2_;


protected:

    // Protected Member Functions

        //- Returns local patchi if patch is in patchIds_ list
        label applyToPatch(const label globalPatchi) const;

        //- Write post-processing info
        virtual void write();


public:

    //- Runtime type information
    TypeName("Bitter");


    // Constructors

        //- Construct from dictionary
        Bitter
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        Bitter(const Bitter<CloudType>& pe);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new Bitter<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~Bitter();


    // Member Functions

        // Evaluation

            //- Pre-evolve hook
            virtual void preEvolve();

            //- Post-patch hook
            virtual void postPatch
            (
                const parcelType& p,
                const polyPatch& pp,
                bool& keepParticle
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Bitter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
