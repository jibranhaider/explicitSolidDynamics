// Acoustic Riemann solver
tC  =
    0.5*(t_M+t_P)
  + (0.5*((fvc::interpolate(Up)*nCn) + (fvc::interpolate(Us)*iMnCn)) & (lm_P - lm_M));

lmC =
    0.5*(lm_M+lm_P)
  + (0.5*((nCn/fvc::interpolate(Up)) + (iMnCn/fvc::interpolate(Us))) & (t_P - t_M));


// Compute boundary values
lm_b.correctBoundaryConditions();
t_b.correctBoundaryConditions();


forAll(mesh.boundary(), patchi)
{
    // Riemann solver for inter-processor boundaries
    if (mesh.boundary()[patchi].coupled())
    {
        vectorField lm_nei = lm.boundaryField()[patchi].patchNeighbourField();

        tensorField P_nei = P.boundaryField()[patchi].patchNeighbourField();
        vectorField Px_nei = Px.boundaryField()[patchi].patchNeighbourField();
        vectorField Py_nei = Py.boundaryField()[patchi].patchNeighbourField();
        vectorField Pz_nei = Pz.boundaryField()[patchi].patchNeighbourField();

        tensorField lmGrad_nei =
            lmGrad.boundaryField()[patchi].patchNeighbourField();

        tensorField PxGrad_nei =
            PxGrad.boundaryField()[patchi].patchNeighbourField();

        tensorField PyGrad_nei =
            PyGrad.boundaryField()[patchi].patchNeighbourField();

        tensorField PzGrad_nei =
            PzGrad.boundaryField()[patchi].patchNeighbourField();

        const vectorField C_nei = C.boundaryField()[patchi].patchNeighbourField();
        const scalarField Up_nei = Up.boundaryField()[patchi].patchNeighbourField();
        const scalarField Us_nei = Us.boundaryField()[patchi].patchNeighbourField();

        forAll (mesh.boundary()[patchi],facei)
        {
            const label& bCellID = mesh.boundaryMesh()[patchi].faceCells()[facei];
            const vector& Cf_ = mesh.Cf().boundaryField()[patchi][facei];

            const vector& lm_M_ =
                lm[bCellID] + (lmGrad[bCellID] & (Cf_ - C[bCellID]));

            const vector& lm_P_ =
                lm_nei[facei] + (lmGrad_nei[facei] & (Cf_ - C_nei[facei]));

            const vector& Px_M_ =
                Px[bCellID] + (PxGrad[bCellID] & (Cf_ - C[bCellID]));

            const vector& Px_P_ =
                    Px_nei[facei] + (PxGrad_nei[facei] & (Cf_ - C_nei[facei]));

            const vector& Py_M_ =
                Py[bCellID] + (PyGrad[bCellID] & (Cf_ - C[bCellID]));

            const vector& Py_P_ =
                Py_nei[facei] + (PyGrad_nei[facei] & (Cf_ - C_nei[facei]));

            const vector& Pz_M_ =
                Pz[bCellID] + (PzGrad[bCellID] & (Cf_ - C[bCellID]));

            const vector& Pz_P_ =
                Pz_nei[facei] + (PzGrad_nei[facei] & (Cf_ - C_nei[facei]));

            const tensor& P_M_ = tensor(Px_M_, Py_M_, Pz_M_);
            const tensor& P_P_ = tensor(Px_P_, Py_P_, Pz_P_);

            const scalar& Up_ = (Up[bCellID] + Up_nei[facei])/2.0;
            const scalar& Us_ = (Us[bCellID] + Us_nei[facei])/2.0;

            const vector& N_ = N.boundaryField()[patchi][facei];
            const tensor& nCn_ = nCn.boundaryField()[patchi][facei];
            const tensor& iMnCn_ = iMnCn.boundaryField()[patchi][facei];

            tC.boundaryFieldRef()[patchi][facei] =
                0.5*((P_M_ + P_P_) & N_)
              + (0.5*((Up_*nCn_) + (Us_*iMnCn_)) & (lm_P_ - lm_M_));

            lmC.boundaryFieldRef()[patchi][facei] =
                0.5*(lm_M_ + lm_P_)
              + 0.5*(((nCn_/Up_) + (iMnCn_/Us_)) & ((P_P_ - P_M_) & N_));
        }
    }

    // Apply boundary conditions
    else
    {
        forAll(mesh.boundary()[patchi], facei)
        {
            lmC.boundaryFieldRef()[patchi][facei] =
                lm_b.boundaryField()[patchi][facei];

            tC.boundaryFieldRef()[patchi][facei] =
                t_b.boundaryField()[patchi][facei];
        }
    }
}