# Finding the phase boundaries of phase diagrams

The overall goal of this project is to find the phase boundaries of coexistence regions within phase diagrams using experimental data. For this specific project presented, the goal of was to find the phase boundary compositions for a two-phase coexistence region varying in composition within the phase diagram of lipid bilayers using experimentally obtained ESR spectra.

## Thermodynamic background

If a multi-component chemical composition is unstable, it can achieve stablility by separating into 2 or more phases with different compositions. Gibbs' phase rule, F = C - P + N, gives the number of intensive thermodynamic variables (i.e. size-invariant) that need to be independently fixed (i.e. thermodynamic degrees of freedom) to uniquely specify the equilibrium state of an open system if the number of chemical components and phases are known, where where P is the number of phases, F is the number of degrees of freedom, C is the number of chemical components, and N is the number of non-compositional variables (like temperature and pressure). In this project, P = 2, C = 2, and N = 0 (because temperature and pressure are constant), which give F = 0. Topologically, the degrees of freedom specify the geometric representation or dimension of the thermodynamic system in the space of the intensive variables.  For example, f = 0 (invariant) is a point, f = 1 (univariant) is a line or curve, and f = 2 (divariant) is a field, plane, or surface. Another intensive variable used is mole fraction, which is the ratio of the amount of a single chemical component over the total amount of all chemical components. A constraint to the system is that the mole fractions sum to unity. In this project, a two-component system is a line of compositions containing a two-phase coexistence region, within which the equilibrium state of the system is uniquely specified as two separated phases each with a different composition. The compositions of each phase are called boundary compositions because they are fixed at each end of the line segment specifying the two-phase region. However, although the boundary compositions do not change within the two-phase region as the total composition changes, the fraction of each phase does change, and these fractions sum to unity.

## Analysis of experimentally obtained ESR spectra

An electron spin resonance (ESR) spectrum is the absorbance of microwave radiation when varing a magnetic field by a free electron spin (i.e. free radical), which is usually attached to a molecule called a spin-probe. An ESR spectrum is independent of phase domain size and is sensitive to phase transitions from one phase to another.  However, a spectrum cannot distinguish individual phases present in coexisting phases without sophisticated spectral analysis. Therefore, the linear superposition model for physical properties of multi-phase systems was applied to analyze the ESR spectra. In general, the linear superposition model of any measurement within a two-phase coexistence region is the weighted sum of the measurements from the two phase boundary compositions, where the weights are 1) the fractions of the measurement from each phase and 2) sum to unity. For application to ESR spectra, this model assumes that the spin-probe partitions between the distinct (possibly submicroscopic or nanoscopic) coexisting phases, and only an insignificant fraction of the probe is at the interface between the phases.  Therefore, because the probe is reporting on the internal physical properties of a phase, the spectrum from a multi-phase system is the linear combination of the spectra from each coexisting phase weighted by the fraction of total spin-probe in that phase.  The process of locating phase boundaries essentially involves finding the optimal basis set of this linear combination with appropriate constraints (such as normalization and non-negativity of the experimental spectra in the absorption mode).  Then, with the proper basis set, the compositions of the coexisting phases can be located; however, so far, this approach can only be applied to binary systems with constant or variable temperature, where the directions (i.e. slopes) of the required thermodynamic tie-lines are known.  Furthermore, the best results are obtained when the spectra from the coexisting phases are significantly different from each other and the spin-probe’s partition coefficient is not much different from unity. 

## Determination of Phase Boundaries in Binary Systems by ESR

Assuming the linear superposition model applies, the method to determine phase boundaries in binary systems from ESR spectra involves the separation of multi-component ESR spectra into a number of linearly independent spectral components.  In the literature this process of component separation is called either self-modeling curve resolution (SMCR) or multivariate curve resolution (MCR).  The method has been applied to uv-vis absorbance spectra of chemical mixtures to determine the individual component spectra; however, it has not been applied to determine compositional phase boundaries in model membrane systems.  Therefore, this project is the application of a version of MCR to ESR spectra along the SPM/DOPC and SPM/Chol binary axes in the SPM/DOPC/Chol lipid system with encouraging results that warrants further development of the method.  

Before analysis, the derivative ESR spectra were integrated into absorbance spectra and then normalized to unit area.   Each absorbance spectrum is converted into an absorbance data vector, which is a discretization of the absorbance intensity versus magnetic field.  The ESR data vectors from sample compositions along the binary axis were organized into a spectral data matrix, S, with the number of rows (M) equal to the number of discrete magnetic field values and the number of columns (N) equal to the number of compositional points sorted by mole fraction from 0 to 1.  This data matrix has the following structure,

 <img width="200" height="40" alt="image" src="https://github.com/user-attachments/assets/10cb6433-57a0-4db9-b56b-57c5107df103" /> (Eqn. 2.1)			
 
where Sα is the matrix of spectral data within the one-phase region of the alpha phase, Sβ is the matrix of spectral data within the one-phase region of the beta phase, and Sγ is the matrix of spectral data within the alpha + beta coexistence region.  A priori the phase boundary spectra that partition the data matrix into these three regions are unknown; however, ordering the columns of the data matrix by increasing mole fraction and knowing that a two-phase region is flanked by one-phase regions determines the matrix structure.  In addition, the constraints imposed on the S matrix were non-negativity (i.e. positive definite) and closure (i.e. normalization),

<img width="400" height="40" alt="image" src="https://github.com/user-attachments/assets/eb945190-8d8a-4876-8c15-c25bca6f2c99" />

<img width="200" height="100" alt="image" src="https://github.com/user-attachments/assets/832ddffb-ebbd-41eb-927d-eb7b16f6729a" /> (Eqns. 2.2)

The non-negativity constraint reflects the fact that absorbance intensities are defined as non-negative scalars; therefore, negative intensities that occur within the noise of the baseline are set to zero.  The closure constraint is just the discrete version of the integral normalization of the absorbance spectra to unit area rescaled by a constant magnetic field interval,  

<img width="500" height="100" alt="image" src="https://github.com/user-attachments/assets/c5205fd0-d47b-41fa-8199-fff29576f6cf" />

Under the linear superposition model, the coexistence spectra, Sγ, are a convex linear combination of the phase boundary spectra,

<img width="300" height="150" alt="image" src="https://github.com/user-attachments/assets/fb645e98-fc24-4f72-8ec4-4beb3917a862" /> (Eqn. 2.3)

with constraints on the coefficients,

<img width="200" height="200" alt="image" src="https://github.com/user-attachments/assets/a64cc91a-a3b4-4664-8079-cc4e4e45bee8" /> (Eqn. 2.4)
						          			   
where fα is the fraction of total spin-probe in the alpha phase, fβ is the fraction of total spin-probe in the beta phase, Sα is the alpha phase boundary spectrum, and Sβ is the beta phase boundary spectrum.  With the definitions,

<img width="300" height="150" alt="image" src="https://github.com/user-attachments/assets/d5d8bc16-940b-45aa-8253-dd8d857368fd" /> (Eqn. 2.5)
					                        
the coexistence spectra can be written in matrix form,

<img width="200" height="60" alt="image" src="https://github.com/user-attachments/assets/1d9ad863-4e9c-4866-aab5-1908ce3cf329" /> (Eqn. 2.6)
								
where B is a M x 2 matrix and f is an 2 x C matrix (C is the number of coexistence spectra).  An approximation to the S data matrix, S*, which includes the full compositional range (i.e. includes the alpha and beta one-phase regions) can be written,

<img width="200" height="60" alt="image" src="https://github.com/user-attachments/assets/24a7a730-6d2b-49f7-aca3-3d5b99adf967" />

where the f* matrix (2 x N)  is the f matrix extended by N−X rows of [1 0] for the alpha phase or [0 1] for the beta phase,

<img width="200" height="200" alt="image" src="https://github.com/user-attachments/assets/a25b720c-2169-47b0-a2cc-96c58a6d8063" /> (Eqn. 2.7)

The assumption in this approximation is that the one-phase spectra are the same as their respective phase boundary spectra.  Because spectral intensities vary with composition in a one-phase region, the validity of this assumption is obviously questionable; however, the linear superposition model only applies to the data near and within the coexistence region, so a method to determine phase boundaries using data from the entire compositional range must systematically and rationally address the flanking one-phase regions.  Furthermore, there is no requirement that the data must explicitly contain the phase boundary spectra or be very near the phase boundary spectra, but there is a requirement that the data overlap both the one-phase regions and the coexistence region.       

Therefore, within the previous formulation where S is known and B, f* are unknown, we seek the solution to the following minimization problem,

<img width="400" height="80" alt="image" src="https://github.com/user-attachments/assets/cfbcd23b-b902-4ac9-99fc-cdec4a0a6ef3" /> (Eqn. 2.8)		

When B and f* are unconstrained, the formal solution to this problem (in the literature called “total least squares”, TLS) is rank reduction through the singular value decomposition (SVD) of the S data matrix.  The rank of the S* matrix is two.  From the SVD of S,

<img width="300" height="60" alt="image" src="https://github.com/user-attachments/assets/825f7407-b057-4f2e-8c91-c7ea100998aa" /> (Eqn. 2.9)

the columns of the U matrix (M x N) are the orthonormal basis vectors (eigenspectra) for the space of data vectors, the W matrix (N x N) is the diagonal matrix of singular values ordered from highest to lowest, and the rows of the V matrix (N x N) are the orthonormal basis vectors of the dual space to the space of data vectors.  The rank two reduction of the S matrix to the solution S* involves keeping the two highest singular values and setting the rest to zero,

<img width="300" height="150" alt="image" src="https://github.com/user-attachments/assets/e5541c83-a910-4b8d-b8a9-d737648539c8" /> (Eqn. 2.10)

where W2 is the modified singular value matrix W and therefore D2 is a N x N matrix that can be truncated to a 2 x N matrix by removing the N−2 rows of zeros.  To conform to the non-negativity and closure constraints (Eqn. 2.2) on the data vectors, the following N number of standard linear least squares problems are solved,

<img width="400" height="80" alt="image" src="https://github.com/user-attachments/assets/a135ddf9-cb5d-4bfa-b4c2-ebd2b4c7b4ad" /> (Eqn. 2.11)

with constraints,

<img width="1024" height="177" alt="image" src="https://github.com/user-attachments/assets/2c7cef3d-2cc4-42a3-a946-1511445b92fa" />

where each Dj2 is a 2 x 1 vector, U1 is the orthonormal basis vector (eigenspectrum) paired with the highest singular value (i.e. first column of U), and U2 is the basis vector paired with the second highest singular value (i.e. second column of U).  The D2 matrix contains the coefficients of the linear combination of the eigenspectra.  The coefficients are the projections of the data vectors onto the eigenspectra.  According to the assumption that the one-phase spectra are the same as the phase boundary spectra, the one-phase spectra should have the same (constant) projected coefficients as the phase boundary spectra.  Also, each data vector has a pair of coefficients, but, because of the normalization constraint, they are dependent, so only the coefficient corresponding to the first eigenspectrum was analyzed for the phase boundary compositions.  In this analysis, the phase boundaries of a two-phase coexistence region along a binary axis can be obtained from a plot of the coefficients of the first eigenspectrum for all N data vectors versus the mole fraction for each data vector.  Visually, the phase boundaries are the points that mark the transition between a flat region and an increasing/decreasing region.  A data-fitting procedure, which employs a simple model relating the D2 and f* matrices, has been developed; however, the method has not been applied in this study because the verification and optimization of the method requires more work.  In addition, the solution to Eqn. 2.8 under the constraints on B and f* (Eqn. 2.2 and Eqn. 2.4, respectively) can be obtained by the method of MCR with alternating least squares, MCR-ALS, but this algorithm has not yet been implemented.   

The two main obstacles to the future development of this method are the systematic, rational modeling of the one-phase regions and the extension to phase diagrams of more the two components.  Regarding the one-phase regions, the relevant question is how to determine if a possible phase boundary spectrum is the true phase boundary spectrum.  This determination must involve both a way to distinguish spectra within the same phase but different composition and to detect a phase transition through changing spectral intensities.  One degenerate case is when spectra are identical within error (i.e. noise) and each, theoretically, are equally likely to be the true phase boundary spectrum.  Another case is when spectra are different but, when looking at a set of data spanning the proposed phase boundary, the changes in the spectral intensities with varying composition within the one-phase region are indistinguishable from the changes in the spectral intensities within the coexistence region caused by the changing fractions of total spin-probe (f α and fβ in Eqn. 2.3).  The fractions of total spin-probe are dependent on the probe’s partition coefficient and the phase’s mole fraction given by the lever rule; therefore, spectral changes through the coexistence region can be modeled because the phase boundary spectra do not change.  However, a suitable model for the changes in spectral intensities within a one-phase region must not merely reproduce the predictions of the model for the coexistence region.  
	
Regarding the extension of the method to phase diagrams of three or more components, its most applicable use would be determining the vertices (i.e. invariant or triple points) of the three-phase triangle.  In this case, the solution would involve the rank three reduction of the data matrix through its SVD.  However, formation of the appropriate data matrix (as in Eqn. 2.1) is nontrivial for a two-dimensional composition space where a compositional point is represented by a vector of mole fractions instead of a scalar mole fraction as for the one-dimensional composition space in a binary system.  In addition, three-phase triangles are not only flanked by three one-phase regions that need to be modeled, but also by three two-phase regions.

## Workflow 
The following workflow assumes that the data matrix S above has been formed from the experimentally obtain ESR spectra. This involves concatenating individual spectra into a matrix, aligning by magnetic field, turning negative values on the baseline to positive (i.e. using absolute value), and then normalization. 

### Load ESR spectra into memory

load()

### Run script

phase_boundary_fit()

