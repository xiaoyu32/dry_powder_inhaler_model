/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "DsmcParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::DsmcParcel<ParcelType>::DsmcParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    U_(vector::zero),
    Ef_(vector::zero),    
    Ei_(0.0),
    typeId_(-1),
    charge_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> U_;
	    is >> Ef_;	    
            Ei_ = readScalar(is);
            typeId_ = readLabel(is);
	    charge_ = readScalar(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&U_),
                sizeof(U_)
              + sizeof(Ef_)
              + sizeof(Ei_)
              + sizeof(typeId_)
	      + sizeof(charge_)	      
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "DsmcParcel<ParcelType>::DsmcParcel"
        "(const Cloud<ParcelType>&, Istream&, bool)"
    );
}


template<class ParcelType>
void Foam::DsmcParcel<ParcelType>::readFields(Cloud<DsmcParcel<ParcelType> >& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<vector> Ef(c.fieldIOobject("Ef", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Ef);
    
    IOField<scalar> Ei(c.fieldIOobject("Ei", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Ei);

    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::MUST_READ));
    c.checkFieldIOobject(c, typeId);

    IOField<scalar> charge(c.fieldIOobject("charge", IOobject::MUST_READ));
    c.checkFieldIOobject(c, charge);

    label i = 0;
    forAllIter(typename Cloud<DsmcParcel<ParcelType> >, c, iter)
    {
        DsmcParcel<ParcelType>& p = iter();

        p.U_ = U[i];
	p.Ef_ = Ef[i];	
        p.Ei_ = Ei[i];
        p.typeId_ = typeId[i];
	p.charge_ = charge[i];
        i++;
    }
}


template<class ParcelType>
void Foam::DsmcParcel<ParcelType>::writeFields
(
    const Cloud<DsmcParcel<ParcelType> >& c
)
{
    ParcelType::writeFields(c);

    label np =  c.size();

    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<vector> Ef(c.fieldIOobject("Ef", IOobject::NO_READ), np);
    IOField<scalar> Ei(c.fieldIOobject("Ei", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<scalar> charge(c.fieldIOobject("charge", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename Cloud<DsmcParcel<ParcelType> >, c, iter)
    {
        const DsmcParcel<ParcelType>& p = iter();

        U[i] = p.U();
        Ef[i] = p.Ef();
        Ei[i] = p.Ei();
        typeId[i] = p.typeId();
        charge[i] = p.charge();
	i++;
    }

    U.write();
    Ef.write();
    Ei.write();
    typeId.write();
    charge.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DsmcParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType& >(p)
            << token::SPACE << p.U()
            << token::SPACE << p.Ef()
            << token::SPACE << p.Ei()
            << token::SPACE << p.typeId()
	    << token::SPACE << p.Ei();
    }
    else
    {
        os  << static_cast<const ParcelType& >(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.U_),
            sizeof(p.U())
          + sizeof(p.Ef())
          + sizeof(p.Ei())
          + sizeof(p.typeId())
	  + sizeof(p.Ei())
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const DsmcParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
