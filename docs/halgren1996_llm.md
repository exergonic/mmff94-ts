# PDF-to-LLM Extraction
**Source:** `halgren1996.pdf`
**Extracted:** 2026-07-24 15:47
**Method:** PyMuPDF text extraction + equation repair

---

### Merck Molecular Force Field. I. Basis,

Form, Scope, Parameterization, and Performance of MMFF94*

## THOMAS A. HALGREN

Department of Molecular Design and Diversity, Merck Research Laboratories, Rahway, New Jersey 07065 Received 20 March, 1995; accepted 31 August, 1995

## ABSTRACT

This article introduces MMFF94, the initial published version of the Merck molecular force field (MMFF). It describes the objectives set for MMFF, the form it takes, and the range of systems to which it applies. This study also outlines the methodology employed in parameterizing MMFF94 and summarizes its performance in reproducing computational and experimental data. Though similar to MM3 in some respects, MMFF94 differs in ways intended to facilitate application to condensed-phase processes in molecular-dynamics simulations. Indeed, MMFF94 seeks to achieve MM3-like accuracy for small molecules in a combined "organic/protein" force field that is equally applicable to proteins and other systems of biological significance. Å second distinguishing feature is that the core portion of MMFF94 has primarily been derived from high-quality computational data-ca. 500° molecular structures optimized at the HF/6-31G* level, 475 structures optimized at the MP2/6-31G* level, 380 MP2/6-31G* structures evaluated at a defined approximation to the MP4SDQ/TZP level, and 1450 structures partly derived from MP2/6-31G* geometries and evaluated at the MP2/TZP level. Å third distinguishing feature is that MMFF94 has been parameterized for a wide variety of chemical systems of interest to organic and medicinal chemists, including many that feature frequently occurring combinations of functional groups for which little, if any, useful experimental data are available. The methodology used in parameterizing MMFF94 represents a fourth distinguishing feature. Rather than using the common "functional group" approach, nearly all MMFF parameters have been determined in a mutually consistent fashion from the full set of available computational data. MMFF94 reproduces the computational data used in its parameterizationo very
well. In addition, MMFF94 reproduces experimental bond lengths (0.014 Å r_{oot} mean square [r_{rns}]), bond angles (1.2°  r_{ms}), vibrational frequencies (61 cm-' *This article includes Supplementary Material available from the author upon request or via the Internet at ftp.wiley.com/public/journals/jcc/supprnat/l7/490 or http://journals.wiley.com/jcc Journal of Computational Chemistry, Vol. 17, Nos. 5&6, 490-519 (1996) 0 1996 by John Wiley & Sons, Inc. CCC 01 92-8651 I96 1050490-30

## MERCK MOLECULAR FORCE FIELD. I

r_{ms}), conformational energies (0.38 kcal/mol r_{ms}), and rotational barriers
(0.39 kcal/mol r_{ms}) very nearly as well as does MM3 for comparable systems. MMFF94 also describes intermolecular interactions in hydrogen-bonded systems in a way that closely parallels that given by the highly regarded OPLS force field. 0 1996 by John Wiley & Sons, Inc. Introduction olecular-mechanics force fields are a crucial M component in the armamentarium used by computational and medicinal chemists for what has become known as “rational drug design.” Early forms of such force fields go back to the work of Hendrickson’ in the 1960s. Many would find par- ticularly noteworthy the work of Allinger and coworkers in developing MMl? MM2: and MM34; that of Kollman and coworkers in developing AMBER5; that of Jorgensen and coworkers in de- veloping OPLS6; that of Karplus and coworkers in developing CHARMM7; and that of Lifson and coworkers in developing CVFF.’ Recent develop- ments, some of which chart important new direc- tions, include the extended CHARMM force field of Momany and Roneg; the DREIDING force field of Mayo et a1.l’; the UFF force field of Rapp&, et al.”; the YETI force field of Vedani and Huhta”; the SHAPES and VALBOND force fields of Landis and  coworker^'^^ 14; the CFF93 force field devel- oped for the Biosym Consortium on Potential En- ergy Functions by Hagler and co~orkersl~,’~; and the MM4 force field of Allinger et al. I7 Like all contemporary force fields, each of the above em- ploys significant physical approximations that limit its accuracy. Moreover, each applies to a different portion of organic/bio-organic chemistry and is derived in a distinctive fashion from a specific selection of data. Given the severity of the approxi- mations that have had to be made, force-field development has been as much an art as a science. One result is that little consensus has been forged either as to what form the force field should take or as to how it should be derived and tested. One can imagine a different situation-one in which essentially all physically significant effects are incorporated accurately into the force field and in which alternative approximations for specific physical terms can be rigorously tested and vali- dated. What makes such a situation imaginable is the steadily increasing computational power avail- able to computational chemists. Such computa- tional power simultaneously makes it possible to employ more complex and more accurate force fields in molecular simulations and to obtain high-quality computational data against which to determine the form of the force field and on which to base its parameterization. Computational theory has already reached the point at which practical ab ab initio methods routinely give results for small- molecule properties that approach experiment in accuracy ’’ while avoiding the large errors that experiment sometimes incurs.lg Moreover, the req- uisite computational data can be obtained rela- tively easily for essentially any system of interest, including many for which no pertinent experimen- tal data are, or are likely to become, available. Δrguably, then, a complex, broadly parameterized force field even now can best be derived from computational data. Efforts to develop improved force fields using computational data and computationally derived insights are already underway in various laborato- ries. In particular, much has been learned about how to model molecular charge distributions accu- rately ’’ and to incorporate induced-dipole effects arising from molecular polarizability.21 These elec- trostatic terms critically affect nonbonded interac- tions. Research based on the use of computational data obtained from ab initio calculations has also been undertaken to better define the form and improve the parameterization of the valence-coor- dinate terms that depend on bond, angle, and torsional distortions.16,22 Particularly noteworthy in the latter regard are the novel fits of the empiri- cal potential energy expression to ab initio relative energies and first and second derivatives em- ployed by Hagler and coworkers in their deriva- tion of CFF93.16 In this series of articles, we report the results of our own initial effort to employ computationally derived information to develop an improved molecular mechanics/dynamics force field. We should note, however, that we have not relied exclusively on computational data. In particular, we have supplemented and extended the range of the core, computationally derived force field, which itself is quite broad, by also parameterizing the force field against a large number of crystallo- graphically determined structures. This combined

## JOURNAL OF COMPUTATIONAL CHEMISTRY

491

## HALGREN

effort has led to what we call the Merck Molecular Force Field (MMFF). We call this initial published version "MMFF94." We should note at the outset that MMFF94 still makes significant approxima- tions in its treatment of important physical interac- tions. Even this version, however, employs compu- tational data of higher quality and broader range than we believe has been utilized in previous efforts. This effort also embodies a particular point of view on what a force field intended for use in bio-organic and pharmaceutical applications should do and on how it should be derived and validated. We expect that growing computational power will soon allow a computationally based approach to be implemented in an even more comprehensive fashion to develop a physically su- perior force field. In the meantime, we believe the performance and range of applicability of MMFF94 warrant its description and use in computational simulations. To this end, we have deposited the parameters as supplementary material in com- puter-readable form.23 Part or all of each parame- ter file is listed in this or in one of the other articles in this ~ e r i e s. ~ ~ — ~ ~ Moreover, we have collaborated with others to implement MMFF93 in CHARMmZ8 and are working to make MMFF94 available in CHARMM,7 the academic version. In addition, MMFF94 is currently being implemented in the BatchMin module of the MacroModel program suite.29 We also hope to be able to distribute OPTI- MOL?' the host molecular-mechanics program for which MMFF94 was developed, through the Quantum Chemistry Program Exchange?* In the next section, we first state the philosophy that underlies the development of MMFF94. We define the form of MMFF94. The fourth section briefly compares the forms of the MMFF94, MM2X, MM2, MM3, and CFF93 force fields. We then de- fine the range of chemical structures for which the computationally derived "core" portion of MMFF has been parameterized and characterize the com- putational data used. Next, we outline the methodology employed in deriving the force field. We then summarize how MMFF94 performs against computational and experimental data in meeting various structural and energetic tests, and subsequently we describe some elements of its implementation in OPTIMOL, CHARMm, and BatchMin. Finally, we summarize this work and sketch some future directions we believe force-field development will take. Subsequent articles in this series will complete the description of MMFF94 by more fully defining: (a) the parameterization of the van der Waals (vdW) and electrostatic representation (part I124); (b) the parameterization of the valence-coordinate terms that determine molecular geometries and vibrational frequencies (part IIIZ); (c) the parame- terization of the torsion terms that then determine conformational energies and torsional barriers (part IV26); and (d) the further extension of MMFF using a combination of experimetnal data extracted from the Cambridge Crystallographic Database, addi- tional computational data, and carefully calibrated empirical rules (part VZ7). Each of these reports also further characterizes the performance of the new force field in reproducing computational and experimental data. In this introductory article, we summarize MMFF94's performance and address the issues that unify its derivation. One further clarification needs to be made. This version of MMFF is primarily intended for use in molecular-dynamics simulations rather than in en- ergy-minimization studies. As a practical matter, the principal distinction between these applica- tions concerns MMFF94's treatment of low-energy inversion barriers at resonance-delocalized tricoor- dinate nitrogen in amides and in such unsaturated amines as vinylamines, anilines, guanines, and nu- cleic-acid bases. In particular, MMFF94 usually gives nonplanar energy-minimized geometries at nitrogen, even for amides, thereby emulating the nonplanar MP2/6-31G*-optimized geometries used in its parameterization. Yet experimen- tal structures, particularly those determined via crystallographic techniques, tend to show planar or nearly planar geometries that reflect time- averaged atomic positions. When used in molecu- lar-dynamics simulations, MMFF94 produces rela- tively flat dynamically averaged structures for such species. Many current pharmaceutical applications, however, rely on energy-minimization methods because of limitations in software and computa- tional resources. For use in such studies, we are developing and intend to soon describe a modified version, currently called "MMFF94s," that yields nearly planar energy-minimized geometries for de- localized trigonal nitrogen?2 The two force fields share most parameters and yield similar, often identical, results for other systems. Basis and Motivation for Formulation of MMFF94 A molecular mechanics/dynamics force field may reasonably be asked to reproduce accurately any or all of a number of molecular properties, 492

## VOL. 17, NOS. 5 &6

## MERCK MOLECULAR FORCE FIELD. I

including the following: molecular geometries. conformational and stereoisomeric energies. torsional barriers and torsion-deformation intermolecular-interaction energies. intermolecular-interaction geometries. vibrational frequencies. heats of formation. energies. Ideally, a single force field would be capable of reproducing these and other molecular properties accurately both in gas-phase and in condensed- phase simulations. Because of their relatively sim- ple construction, however, current force fields nec- essarily make a variety of compromises. Here we discuss the choices we have made in developing

## MMFF94.

A pivotal application for MMFF94, from which a number of constraints on its design and imple- mentation follow, is the study of receptor-ligand interactions involving proteins or nucleic acids as receptors and a wide range of chemical structures as ligands. For quantitative study, the force field must be able to describe the ligand and receptor properly in isolation as well as when bound. For these purposes, molecular geometries need to be good, but conformational energies are crucial if the force field is to avoid modeling the wrong con- former of the ligand (or receptor) upon binding or giving an erroneous estimate of the energetic cost of adopting the detailed conformation required for binding. To assess these aspects properly, the force field must be able to locate conformational minima accurately and describe intervening torsional pro- files and barriers reasonably well. At least equally importandy, intermolecular- interaction energies (and, to a lesser extent, geome- tries) must also be described accurately. In con- trast, vibrational frequencies should be reasonably accurate, but spectroscopic precision is unlikely to be required. Thus, fine details of vibrational spec- tra, such as the splitting of high-frequency modes for bond stretching or angle bending, are unlikely to appreciably affect the differential free energy of binding to a macromolecular receptor of one lig- and relative to another. Finally, though heats of formation are crucial in some applications, they are not required to understand differences in free energies of binding and are not addressed in

## MMFF94.

To be routinely and reliably useful in pharma- ceutical, bio-organic and chemical applications, MMFF94 would need to be able to handle most organic structural types represented in the Merck Index3 or the Fine Chemicals Directory. M This broad intended range of application places significant requirements on the data to be used in the param- eterization of the force field. We note in this regard that Allinger and coworkers have crafted a series of highly regarded molecular-mechanics force fields based primarily on the meticulous examina- tion and careful selection of good quality experi- mental data,2-4,17 and that other force fields such as AMBER5 and CHARMM7 have also been pa- rameterized mainly against experimental data. This approach, however, could not be used to derive MMFF94, for two reasons. First, the location, selec- tion, and extraction of good experimental data is a highly time-consuming enterprise and requires a degree of expertise we lack. Second, and more importandy, high-quality experimental data, par- ticularly for conformational and intermolecular-in- teraction energies, are unavailable for a great many of the chemical structures MMFF94 must handle. For these reasons, the core portion of MMFF94, on which we focus here, has been derived primar- ily from ab initio data (though experimental data have been liberally employed in its validation). An especially cogent argument for the use of such computational data has recently been offered by Hagler and coworkers.lb" In a novel and notewor- thy departure from previous practice, these work- ers employed data for molecular dipole moments, relative energies, and Cartesian first and second derivatives obtained from HF/6-31G* calculations to characterize the quantum mechanical energy surface used to derive the QMFF (quantum-mech- anical force field) predecessor of CFF93, the Biosym Consortium force field. The approach we have taken in deriving MMF94 is, in part, patterned after theirs. Both, for example, employ the power- ful Consortium program PROBE% to derive force constants for terms related to bond stretching and angle bending from the information on the curva- ture of the quantum mechanical surface contained in the HF/6-31G* second derivatives. However, the two approaches also differ in a number of ways that may materially affect their performance in molecular sim~dations."-~~ The derivation of a force field from computa- tional data would be straightforward if we wished to describe only gas-phase systems. However,

## JOURNAL OF COMPUTATIONAL CHEMISTRY

493

## HALGREN

while many, and perhaps most, of the processes we wish to model occur in condensed phases, MMFF94 accounts for the effects of molecular po- larizability only in a limited way. These effects, for example, cause the dipole moment of water to rise from a gas-phase value of 1.85 D to a mean value of -2.4 D in aqueous s ~ l u t i o n. ~ ' ~ ~ ~ ~ Clearly, a condensed-phase simulation that uses a gas-phase dipole moment for water would seriously under- estimate electrostatic interactions and would be expected to yield poor computational proper- ties.37 Consequently, MMFF94, like OPLS6 and other current force fields intended for use in con- densed-phase simulations, employs effect pair potentials37 that reflect, in an averaged sense, the enhancement of the charge distribution due to molecular polarizability. Especially careful attention must be given to the partitioning between electrostatic, van der Waals (vdW), and torsional  interaction^.^' Our approach begins by choosing the vdW representation as pre- viously defined39 and the electrostatic representa- tion from fits to scaled (enhanced by HF/6-31G* 40 molecular dipole moments. To prop- erly describe hydrogen-bonding interactions, we then adjust key vdW and electrostatic parameters to better fit scaled intermolecular-interaction ener- gies and geometries obtained from HF/6-31G* calculation^.^^ Last, we derive the torsion terms to fit the ab initio gas-phase conformational data. Conveniently, the HF/6-31G* level of theory con- sistently overestimates gas-phase dipole moments for organic compounds.6b'1xa For water, it gives a calculated dipole moment of 2.20 D,42 or 2.42 D after 10% enhancement, close to the previously cited mean value of -2.4 D found in aqueous solution. This use of scaled HF/6-31G* interac- tion energies and geometries allows MMFF to be parameterized in a straightforward manner that seeks to ensure that a proper balance between solvent-solvent, solvent-solute, and solute-solute interactions is achieved. The quality of this balance is crucial for accurately describing aqueous solva- tion and the energetics of host/guest binding in aqueous solution. We have also explored the use of higher level ab initio  calculation^,^^ but have not found an alternative approach that appears prefer- able. As noted in the Introduction, one further choice we made in developing MMFF94 was to derive a force field explicitly intended for use in molecular-dynamics simulations. Å modified ver- sion more suitable for use in energy minimization studies (MMFF94s) is also being developed.~' Form of the Merck Molecular Force Field The MMFF94 energy expression can be written as:

$$
E = ΣEB_{ij}
+ ΣEA_{ijk} + ΣEBA_{ijk}
+ ΣEOOP_{ijkl} + ΣET_{ijkl}
+ ΣEvdW_{ij} + ΣEQ_{ij}
(1) where the seven constituent terms are defined as shown below. In each case, the cited numerical constant is such that the deformation or interaction energy is expressed in kilocalories per mole when distances and angles are measured in angstroms and in degrees, respectively. In the notation that follows, we adopt the con- vention that a specific atom involved in a force- field interaction is designated by i, j, k,... and that the corresponding numeric MMFF atom type is designated by I, I, K,... This notation makes explicit, for example, that the force constant k b,, and the reference bond length r₀_{IJ} for the i-j bond in eq. (2) depend on the associated MMFF atom
types 1 and J, whereas the bond distance, r_{ij}, depends on the atomic coordinates. In certain in- stances, the parameters depend only on the atomic species for atoms i, j, k,...; in such cases, we still use capital letters, but explicitly note the actual dependence in the text.
$$


## BOND STRETCHING

MMFF94 employs the quartic function:

$$
EB_{ij} = 143.9325-Δr² k_{b,IJ}/2
x(1 + cs Δr_{ij} + 7/12cs² Δr_{ij}²)
$$

(2)
where kb_{IJ} is the force coptant (md/&, Δr_{ij} =
r_{ij} — rp, is the difference (Å) between $ctual and reference bond lengths, and cs = -2 A-' is the "cubic-stretch' constant. This function corre- sponds to an expansion through fourth oorder of a Morse function with an "alpha" of 2 A-1.43 Re- sults published in a recent high-level ab initio study44 show this value for alpha to be a represen- tative one. Special sets of reference bond lengths 494

## VOL. 17, NOS. 5 & 6

## MERCK MOLECULAR FORCE FIELD. I

and force constants are employed for "conjugated single bonds," such as those found in butadiene and biphenyl, as well as for certain other single bonds between sp- or spz-hybridized atoms.z5

## ANGLE BENDING

MMFF94 normally uses the cubic expansion:

$$
EA_{ijk} = 0.043844 ka_{IJK}/2 Δθ_{ijk}(1 + cb Δθ_{ijk}) (3) 2 where ka_{IJK}, is the force constant (md A/radz),
$$

= θ_{ijk} — θ₀_{IJK} is the difference between ac- tual and reference bond angles (degrees), and cb = -0.007 deg-' (or, more precisely, -0.4 rad-'1 is the "cubic-bend" constant. Special sets of parame- ters are used for angles that involve delocalized single bonds and/or occur in small ringsE For linear or near-linear bond angles, MMFF94 em- ploys the well-behaved form used in DREIDING" and UFF":

$$
EAjjk = 143.9325ka_{IJ}(l + cos θ_{ijk})
$$

(4)

## STRETCH-BEND INTERACTIONS

MMFF94 employs the form:

$$
EBAijk = 2.51210(kba_{IJK} Δr_{ij} + kba_{KJI} Δ_{rkj}) Δθ_{ijk}
(5) where kbu,,, and kba_{IJK} are force constants (md/rad) that couple the z-j and k-j stretches to the i-j-k bend, and Å T and Å 8 are as defined above. Stretch-bend interactions are omitted when eq. (4) is used for bond angles.
$$


## OUT-OF-PLANE BENDING AT TRICOORDINATE

## CENTERS

MMFF94 uses the form: where koop_{IJK}: is the force constant (md A/rad2)
and χ_{ijkl} is the Wilson angle45 (degrees) between the bond j-1 and the plane i-j-k. The three angles that arise at a given center, j, are all assigned the
same koop_{IJK}: force constant; the "in-plane" an- gles use "normal" bond angles and are described by eq. (3). For trigonal nonplanar centers, this formulation allows angle-bending reference values that average less than 120°  to be used to make the center pyramidal; the out-of-plane term can then be employed to improve the fit to the inversion barrier.

## TORSION INTERACTIONS

MMFF94 uses the threefold representation em- ployed in MM2 and MM3 and MM3: where Q, is the i-j-k-Z torsion angle:

$$
ET_{ijkl} = 0.5(V1(1 -k cos φ) + v2(1 — cos 2φ) + v3(1 + cos 3φ)) (7)
$$

The constants V_{ij} V_{ij} and V, depend on the atom types I, J, K, and L for atoms i, j, k, and I, where i-j, j-k, and k-1 are bonded pairs. Torsion interac- tions within four-membered rings and saturated five-membered ringsz6 are given special torsion constants, as are interactions in which either the central or a wing bond is a single bond between two atoms that are capable of participating in multiple or aromatic bonds?6 The former situation occurs, for example, in biphenyl, butadiene, and styrene.

## VAN DER WAALS INTERACTIONS

MMFF employs the recently developed "Buffered-14-7°  form39; the terminology derives from the formal 14th and 7th power dependencies for the repulsive and attractive terms that would be obtained if the RT, "buffering constants" in the denominators were deleted. The form of the poten- tial is shown in eq. (8):
(8) This form is used in conjunction with an expres- sion that relates the minimum-energy separation RT, to the atomic polarizability aI [eq. (911, with specially formulated combination rules [eqs. (10) and (ll)], and with a Slater-Kirkwood expression
for the well depth E_{ij} [eq. (1211:

## RT, =

(9)

## JOURNAL OF COMPUTATIONAL CHEMISTRY

495

## HALGREN

As described elsewhere,’~ modified values of Ry, and c1, are used to describe hydrogen-bonding interactions. Van der Waals and electrostatic inter- actions are included whenever atoms i and j be- long to separate domains or are separated by three or more chemical bonds; 1,4-vdW interactions are not differentially scaled in MMFF94.

## ELECTROSTATIC INTERACTIONS

MMFF94 uses the buffered coulombic form:

$$
EQ_{ij} = 332.0716q_i q_j/(D(R_{ij} + 6)”)
$$

(13)
where 9, and q_i are partial atomic charges, R_{ij} is
the internuclear separation in A, 6 = 0.05 Å is the “electrostatic buffering” constant, and D is the “dielectric constant.” Normally, the exponent n is taken as 1, though use of a distance-dependent dielectric constant (n = 2) is also supported. In MMFF94, 1,4-electrostatic interactions are scaled by a factor of 0.75.26 The distance buffering, where 6 > 0, prevents infinite attractive electrostatic en- ergies from overwhelming the finite repulsive vdW interaction contained in eq. (8) as oppositely charged atomic centers coalesce. The partial atomic charges q1 used in eq. (13) are constructed from initial full or fractional ”for- mal atomic charges” qp (usually zero, but, e.g., + 1/3 for guanidinium nitrogen) by adding contri- butions from bond charge increments wKl that describe the polarity of the bonds to atom i from attached atoms k. Specifically, MMFF94 computes 41 as
where w_{KI} = — wIK. The procedure used to assign the qp is specified in part V of our st~dy.’~ Comparison of MMFF94’s Functional Form to MM2, MM2X, MM3, and CFF93 MMFF94 closely resembles MM2 and MM3, as well as MM2X, our previous generation force in functional form. For bond stretching, MMFF94 and MM3 each use a quartic expansion in which the cubic and quartic force constants are related to the quadratic force constants in a prede- termined way. Each thereby avoids the “cubic stretch” catastrophe, in which progressive elonga- tion of a chemical bond eventually drives the MM2 or MM2X energy to negative infinity. This catas- trophe is circumvented in MM2X’s implementa- tion at the cost of additional complexity in the computer code that might prove troublesome in molecular-dynamics applications. MMFF94 and MM3 employ anharmonic angle bending, an in- trinsically better representation than the simple quadratic form used in MM2 and MM2X, though MMFF truncates its representation at the cubic term.47 Moreover, trigonal centers are handled dif- ferently in MMFF94 to allow out-of-plane terms to be used for certain centers that have nonplanar equilibrium ge~metries.’~ Centers having linear idealized bond angles are also handled differently, through eq. (4). The same forms for stretch-bend and torsion interaction are used in all four force fields. MMFF94 currently omits MM3’s bond-tor- sion and bend-bend terms; bond-torsion terms may be included later in certain cases, as may angle-torsion terms. MM3 also includes elec- tronegativity-related adjustments to reference bond lengths that MMFF94 omits but that in certain cases can be ~ignificant.’~,’~ The most important differences between MMFF94 and the MM2, MM2X, and MM3 force fields arise in the description of nonbonded inter- actions. In particular, MMFF94 uses a “buffered” expression for vdW interactions and employs novel combination rules for the vdW parameters. Unlike MM2X>6 MMFF94 properly treats intramolecular and intermolecular electrostatic interactions in the same manner. Moreover, MMFF94 normally uti- lizes a unit dielectric constant, thereby allow- ing the force field to be applied without modifica- tion to condensed-phase simulations employing explicit solvent. In addition, like AMBER,5 CHARMM,7 CVFF?, CFF93,I6 and most other force fields used in molecular-dynamics simulations, MMFF94 describes hydrogen-bonding interactions as being essentially electrostatic in nature, whereas MM2 (1987 parameters and later) and MM3 obtain up to 6 kcal/mol of stabilization energy4 from an attractive Exp-6 term. CFF93 and MMFF94 both use a quartic expan- sion for bond stretching, treat stretch-bend inter- actions in the same way, and employ equivalent representations for torsion interactions. I6 Both, in addition, define partial atomic charges in terms of bond charge increments, describe electrostatic in- teractions solely in terms of charge-charge interac- tions (i.e., avoid special “hydrogen-bond” terms5r7), and use novel vdW combination rules39,49 in conjunction with a vdW potential (Lennard-Jones 9-6 or Buf-14-7) that differs from the more commonly used Lennard-Jones 12-6 496

## VOL. 17, NOS. 5 & 6

## MERCK MOLECULAR FORCE FIELD. I

form. Minor differences include MMFF94's use of a cubic rather than quartic expansion for angle bending and its use of a special form [eq. (4)] to describe "linear bond angles." One major differ- ence is that CFF93 includes many more "cross terms." These terms allow CFF93 to describe cer- tain elements of geometry more accurately, an example being the elongatio? of a conjugated sin- gle bond by as much as 0.1 Å upon bond rotation. The additional cross terms also allow CFF93 to reproduce vibrational spectra more accurately than can MMFF94. As noted earlier, however, MMFF94 does not seek to predict vibrational spectra to high accuracy but rather to describe conformational and intermolecular interaction energies as well as pos- sible. With respect to these latter considerations, major differences in parameterization that may materially affect performance arise from differ- ences in the data and methodology ~ s e d? ~ — ' ~ Survey of Systems for Which MMFF94 is Parameterized To define the range of organic and bio-organic systems covered in the core parameterization of MMFF94, the set of compounds and molecular conformations for which optimized MP2/6-31G* geometries have been employed are listed in Table I. Each structure is identified by a five-char- acter conformational index of the form "XYNMc," where "XY" defines the compound class (e.g., AM for amides and related compounds), "NM" speci- fies the compound number within the class, and "c" identifies the conformation. The conforma- tional designations "u" through "i" correspond to equilibrium conformers; " j" through "2°  indicate conformations optimized while holding a particu- lar torsion angle fixed, except that "t" and some- times "s" usually denote a symmetry-determined conformational transition state. To characterize the conformation, the compound name is followed, where appropriate, by a brief description of the geometry. Among "monofunctional" chemical families, MMFF94 has been parameterized for alkanes, alkenes, alcohols, phenols, ethers, aldehydes, ke- tones, ketals, acetals, hemiketals, hemiacetals, amines, amides, peptide analogs, ureas, imides, carboxylic acids, esters, carboxylate anions, ammo- nium cations, thiols, mercaptans, disulfides, halides (chlorides and fluorides), imines, iminium

## TABLE 1.

Conformers Considered in Parameterizing Core

## M M FF94.'


$$
Amides and peptide analogs AMOla -formamide AMOlt -formamide, N planar AM02a -N-methylformamide, cis AM02b -N-methylformamide, trans AM02j -N-methylformamide, trans, h-c-n-c = 60° AM02k -N-methylformamide, trans, h-c-n-c = 30° AM021 -N-methylformamide, trans, h-c-n-c = 0° h-n-c=o = 115° h-n-c= o = 120° h-n-c= o = 125° h-n-c= o = 130° h-n-c= o = 55° h-n-c= o = 60° h-n-c= o = 65° h-n-c= o = 70° AM02s -N-methylformamide, AM02t -N-methylformamide, AM02u -N-methylformamide, AM02v -N-methylformamide, AM02n-N-methylformamide, AM02x -N-methylformamide, AM02y -N-methylformamide, AM02z -N-methylformamide, AM03a -acetamide AMO3t -acetamide, N planar — anti ts, — anti ts, — anti ts, — anti ts, — syn ts, — syn ts, — syn ts, — syn ts, AM04a -N-methylacetamide, trans AM04b -N-methylactamide, cis AM04j -N-methylacetamide, trans, AM04k -N-methylacetamide, trans, AM041 -N-methylacetamide, trans, AM04m-N-methylacetamide, cis, AM04s -N-methylacetamide, cis, AM04t -N-methylacetamide, trans, AM05a -N, N-dimethylformamide AM06a -urea, puckered AM06t -urea, planar AM07a -N-formylformamide, both o= c-n-h cis AM07b -N-formylformamide, o= c-n-h cis, trans AM08a -formylglycinamide AM09a -glycine dipeptide analog, C7 AM09b -glycine dipeptide analog, C5 h-c-c= o = 0° h-c-c= o = 30° h-c-c= o = 60° h-c-c= o = 0° N planar N planar (Continues on next page)
$$


## JOURNAL OF COMPUTATIONAL CHEMISTRY

497

## HALG R EN

## TABLE I.

(continued)

## AMO9s — glycine dipeptide analog, C7, N planar

## AMO9t — glycine dipeptide analog, C5, N planar AMlOa — alanine dipeptide analog, C7eq AMlOb — alanine dipeptide analog, C5 AMlOc — alanine dipeptide analog, C7ax AMlOd — alanine dipeptide analog, a' AMlOe — alanine dipeptide analog, p2 AMlOf — alanine dipeptide analog, aL

## AM1 1 a — propionamide, c — c — c — n anti

## AM1 2a — N-ethylformamide, c- c- n -c gauche

## AM1 2a — N-ethylformamide, c — c — n — c = 180°

## AM1 3a — N-OH, N-methylacetamide, o — n — c= o trans

## AM1 3b — N-OH,N-methylacetamide, o-n -c= o cis

## AM1 3s — ## N-OH NMA, o — n — c= o trans, N planar

## AM1 3t — ## N-OH NMA, o — n — c= o cis, N planar AM14a — ## N-OH,N — Et acetamide, o — n — c= o trans, c-c-n-ogauche o — n — c= o trans, c-c- n — o trans AM14c — ## N-OH,N — Et acetamide, 0- n -c= o cis, c-c-n-c(=o) skew AM14d — ## N-OH, N — Et acetamide, o — n — c= o cis, c — c- n — c(=o) gauche AM15a — N-OH, N-Me propionamide, 0- n -c= o trans, c- c- c= o cis AM15b — N-OH, N-Me propionamide, o — n — c= o trans, c- c- c= o skew AM15c — N-OH, N-Me propionamide, o — n — c= o cis, c-c-c=o cis

## AM1 5d — N-OH, N-Me propionamide, o-n -c= o cis, c- c- c= o skew AM14b — ## N-OH,N — Et acetamide, AM16a — glycine dipeptide, C7 AM16b — glycine dipeptide, C5 AM17a --lanine dipeptide, C7eq AM17b --lanine dipeptide, C5 AM17c --lanine dipeptide, C7ax AM17d --lanine dipeptide, a'

## AM1 7e — alanine dipeptide, p2 AM17f — alanine dipeptide, aL

## TABLE 1.

(continued) Carboxylate anions

## AN01 a — formate anion AN02a -acetate ion AN03a — propionate anion AN04a — propenoate anion Δromatic and heteroaromatic compounds

## ARO1 a — benzene AR02a — pyridine AR03a — pyrimidine AR04a — pyridazine AR05a — 1,3,5-triazine

## ARO6a — pyrrole AR07a -furan AR08a — thiophene

## ARO9a — imidazole ARlOa — pyrazole

## AR11 a — 1,2,4-triazoIe AR12a — 1,2,3,4-tetrazoIe (N1 tautomer) AR13a — 1,2,3,5-tetrazole (N2 tautomer) AR14a — oxazole AR15a — isoxazole AR16a — 1,3,4-0xadiazole AR17a — 1,2,4-0xadiazole ARi 8a — thiazole AR19a — isothiazole AR20a — 1,3,4-thiadiazole

## AR21 a — pyridine N-oxide AR22a -toluene AR23a — ethylbenzene, c — c — c — c skew AR23t — theylbenzene, c — c — c — c = 0° AR24a — N-ethylpyrrole, c — n — ch, — ch, ca. 90°° AR25a — indole Carboxylic acids CAOl a — methanoic acid, o= c — o — h cis CAOl b — methanoic acid, o= c — o — h trans CA02a — ethanoic acid, o= c — o — h cis CA02b — ethanoic acid, o= c — o — h trans CA03a — propanoic acid, c — c — c= o cis

## CAO3b — propanoic acid, c — c — c= o skew CA04a — glyoxalic acid, o= c — c= o trans, o= c — o — h cis

$$
CA04b -glyoxalicacid, o=c-c= o trans, o=c-o-h trans (h-bond) CA05a — glycolic acid, o= c — c — o cis (h-bond)
$$


$$
CA05b -glycolic acid, o= c — c — o skew (h-bond)
$$


## CAO6a — benzoic acid CA07a — propenoic acid, c= c — c= o trans CA07b — propenoic acid, c= c — c= 0 cis CA08a -oxalicacid, o=c-c=o trans, both o= c — o — h trans CA08b -oxalicacid, o=c-c=o trans, both o= c-0- h cis (Continues on next page) 498

## VOL. 17, NOS. 5&6

## MERCK MOLECULAR FORCE FIELD. I

## TABLE 1.

(continued)

## CAO8c — oxalic acid, o= c — c= o trans, one o= c-0- h trans

## CAO9a — pyruvic acid, o= c — c= 0, o = c — o — h trans (h-bond)

## CAO9b — pyruvic acid, o= c — c= 0, o=c-0-hcis Carboxylic acid esters CEOla -methyl formate, o= c — o — c cis

## CEO1 b — methyl formate, o= c — o — c trans

## CEO1 j — methyl formate, o= c — o — c trans,

## CEO1 k — methyl formate, o= c — o — c trans,

## CEO1 I — methyl formate, o= c — o — c trans, CE02a — methyl acetate, o= c — o — c trans CE02b — methyl acetate, o= c — o — c cis CE05a — vinyl formate, o= c — o — c cis, c= c — o — c trans CE05b — vinyl formate, o= c- o — c cis, c= c-0- c cis CEOGa -ethyl formate, o= c — o — c cis, c — o — c — c anti

## CEOGb — ethyt formate, o= c — o — c cis, c-0-c-c gauche CE07a — isopropyl formate, c — c — o — c = g, a CE07b — isopropyl formate, c — c — o — c = g, g CE08a — phenyl acetate CE08b — phenyl acetate, c — c(= 0) — o — c — 90°

## CEOSa — propiolactone CElOa — methyl glycolate, o= c- c — o cis (h-bond)

## CE1 Ob — methyl glycolate, Conjugated systems CJOla — 1,3-butadiene, gauche

## C J O l b — 1,3-butadiene, s-trans CJOl t — 1,3-butadiene, c= c — c= c = 0° CJ02a — 2-methyl-l,3-butadiene, gauche CJ02b — 2-methyl-l,3-butadiene, s-trans CJ03a — 2-methyl-but-1 -ene-3-one, c= c — c= o cis CJ03b — 2-methyl-but-1 -ene3-one, c= c — c= o trans CJ04a — 2-methylpropenamide, c= c — c= o cis

## CJO4b — 2-rnethylpropenamide, c= c — c — o skew CJ05a — propenamide, c= c — c= o cis CJ05b — propenamide, c= c — c= o skew

## CJO6a — but-1 -ene-3-one, c= c — c= o cis

## CJOGb — but-1 -ene-3-one, c= c — c= o trans CJ07a — acrolein, c= c — c= o cis CJ07b — acrolein, c= c — c= o trans CJ08a — 2-methylpropenal, c= c — c= o cis CJ08b — 2-methylpropenal, c= c — c= o trans

## CJO9a — 2-methylpropenoic acid, c= c — c= o trans

## CJO9b — 2-methylpropenoic acid, c= c — c= o cis

## h-C-0-C

= 180°

## h-C-0-C

= 150°

## h-C-0-C

= 120° o = c — c — o skew (h-bond)

## TABLE 1.

(continued)

## CJ1 Oa — acetophenone

## CJ11 a — styrene CJ12a — 2-phenylpropene CJ12j — 2-phenylpropene, framework planar CJ13a — 1,3-pentadiene, s-trans, c — c= c — c trans CJ13b — 1,3-pentadiene, gauche, c — c= c — c trans CJ13c — 1,3-pentadiene, s-trans, c — c= c — c cis CJ14a — cyclopentadiene Aldehydes and ketones COO1 a -formaldehyde C002a -acetaldehyde C003a — propionaldehyde, c — c — c= o cis C003b — propionaldehyde, c — c — c= o C004a -acetone C005a — butanone c — c — c — o = 0° C005b — butanone c — c — c — o skew C005j — butanone, c — c — c= o = 0° C005k — butanone, c — c — c= o = 30°

## COO51 — butanone, c — c — c= o = 60° C006m- butanone, c — c — c= o = 90° C005n — butanone, c — c — c= o = 120°

## COO50 — butanone, c — c — c= o = 150° C005p — butanone, c — c — c= o = 180°

## COO6a — methyl isopropyl ketone, o= c — c(ch,), — h trans C006b — methyl isopropyl ketone, o= c — c(ch,), — h cis C007a — butyraldehyde, c — c — c — c anti C007b — butyraldehyde, c — c — c — c gauche C008a — but-3-enal c= c — c — c skew, c — c — c= o cis C008b — but-3-enal c= c — c — c skew -, C008c — but-3-enal c= c — c — c skew+, C009a — 3-methyl-but-3-ena1, c= c — c — c skew, C009b — 3-methyl-but-3-ena1, c= c — c — c skew, COlOj -isobutyraldehyde, h-c(=o)-c- h = 0° COlOk — isobutyraldehyde, h — c(=o) — c — h = 30° COlOl -isobutyraldehyde, h-c(=o)-c-h = 60°

## CO1 Om-

isobutyraldehyde,

## CO1 On — isobutyraldehyde, C0100 — isobutyraldehyde,

## CO1 Op — isobutyraldehyde, C011 a — cyclobutanone C011 t — cyclobutanone, planar c- c- c= o skew+ c- c- c= o skew+ c- c- c= o cis c- c- c= o skew h-c(=o)-C-h=90° h--C(=o)-C-h = 120° h-C(=O)-C-h = 150° h-C(=O)-C-h = 180° (Continues on next page)

## JOURNAL OF COMPUTATIONAL CHEMISTRY

499

## HALGREN

## TABLE 1.

(continued) C012a — 2-formylpropana1, o — c — c — c(= 0) anti C012b — 2-formylpropanal, C013a -4-oxobutanal, o=c-c-c, C013b -4-oxobutanal, o= c-c-c, C014a — 2,3-butanedione, c — c — c — c trans C014t — 2,3-butanedione, c — c — c — c cis Halides HLOl a — fluoromethane HL02a — difluoromethane HL03a — 1,l -difluorethane HL04a — 1,2-difIuoroethane, f — c — c -f anti HL04b — 1,2-difluoroethane, f — c — c -f gauche HL05a — 1,e-dichloroethane, cI — c — c — cI anti HL05b — 1,Bdichloroethane, cI — c — c — cI gauche

## HLO6a — 1,l,l-trifluoroethane HL07a — l,l,l-trichloroethane HL08a — chlorocyclobutane

## HLO8j — chlorocyclobutane, planar

## HLO9a — fluoropropane, c — c — c — f anti

## HLO9b — fluoropropane, c — c — c — f gauche

## HL1 Oa — chloropropane, c — c — c — cl anti

## HL1 Ob — chloropropane, c — c — c — cI gauche Imines, guanadines, and amidines IMOla — formamidine, h — n = c — n cis,

## IM01 b — formamidine, h — n= c — n anti, lMOl t — formamidine, h — n= c — n cis, N planar IM02a — N-methylformaldehydeimine, h — c — n= c cis IM02t — N-methylformaldehydeimine, h-c-n=c= 180° IM03a — formaldehydeimine IM04a — N-methylformamidine, n — c= n — c cis, N puckered IM04b — N-methylformamidine, n — c= n — c trans, N puckered IM04t — N-methylformamidine, n — c= n — c cis, N planar IM05a — guanidine, N puckered IM05t — guanidine, planar IM06a — N2-methylguanidine, N puckered IM06t — N2-methylguanidine, N planar 1M07a — butadiene Schiff base, c= c — c= n IM07b — butadiene Schiff base, c= c — c= n s-cis, 0-c-c-c(=o) gauche c-c-c=o cis, c-c-c-ctrans c-c-c=o cis, c-c-c-c gauche N puckered N puckered s-trans, h — n = c — c cis h — n= c-c trans 500

## TABLE I.

(continued) Ketals, acetals, and hemiacetals KT02a — 2-methoxytetrahydropyran, equatorial KT02b — 2-methoxytetrahydropyran, axial, KT03a — 2,4 dioxapentane, c — o — c — o g +, KT03b — 2,4 dioxapentane, c — o — c — o g, KT04a — 2,5-dimethyl-l,3-dioxane (5-equatorial) KT04b — 2,5-dimethyl-l,3-dioxane (5-axial) KT05a — methoxymethanol, c — o — c — o g +, KT05b — methoxymethanol, c — o — c — o g +, KT05c — methoxymethanol, c — o — c — o g, Cations NCOl a -ammonium cation NC02a — N-methylamine cation NC03a — N-ethylamine cation NC03t — N-ethylamine cation, h — n — c — c = 0° NC04a — N,N-dimethylamine cation NC05a — N-propylamine cation, c-c-c-n gauche NC05b — N-propylamine cation, c — c — c — n anti

## NCO6a — guanidine cation NC07a — ethylguanidine cation, c — c — n = c anti NC07b — ethylguanidine cation, c — c — n= c gauche NC08a — formamidine cation

## NCO9a — methylguanidine cation

## NC1 Oa — N-methylformaldehydeimine cation NCl 1 a — N-methylformamidine cation, c-n-c=ncis

## NC11 b — N-methylformamidine cation, c — n — c= n trans NC12a — imidazolium cation NC13a — formaldehydeimine cation NC14a — t-butylamine cation

## OCO1 a — hydronium ion Amines

## NHO1 a — methylamine NH02a — propylamine, c — c — c — n anti NH02b — propylamine, c — c — c — n gauche NH03a — isopropylamine, C1 (c- h gauche to n -Ip) NH03b — isopropylamine, Cs (c- h anti to n — Ip)

## NHO3j — isopropylamine, h — c — n — h = 120° NH03k — isopropylamine, h — c — n — h = 150° NH031 — isopropylamine, h — c — n — h = 180° NHO3m- isopropylamine, h — c — n — h = 21 0° NH03n — isopropylamine, h — c — n — h = 240° NH03p — isopropylamine, h — c — n — h = 270° NH03p — isopropylamine, h — c — n — h = 300° me-0-c-canti 0-c-0-cg+ 0-c-0-ca o-c-o-hgi- 0-C-0-hg- o-c-0-ha (Continues on next page)

## VOL. 17, NOS. 5 & 6

## MERCK MOLECULAR FORCE FIELD. I

## TABLE 1.

(continued) NH04a — cyclohexylamine, equatorial NH04b — cyclohexylamine, axial NH05a — dimethylamine

## NHO6a — azetidine, n — h equatorial

## NHO6j — azetidine, ring planar NH07a — piperidine, n — h equatorial NH07b — piperidine, n — h axial NH08a — trimethylamine

## NHO9a — N-methylpiperidine, equatorial

## NHO9b — N-methylpiperidine, axial NHlOa -ammonia NHlOt -ammonia, planar NHll a — ethylamine, c — c — n — Ip gauche NHll b — ethylamine, c — c — n — Ip anti NH12a — f-butylamine NH13a — vinylamine NH14a -aniline, N puckered NH14t -aniline, planar NH15a — pyrrolidine, n — h equatorial NH15j — pyrrolidine, ring planar NH16a — 3-aminopropene, c= c — c — n skew, c-c-n-lpgauche NH16b — 3-aminopropene, c= c- c — n cis, c- c- n — Ip gauche NH16c — 3-aminopropene, c= c — c — n skew, c- c- n — Ip anti NH17a — 2-me,3-aminopropene, c= c — c — n skew, c-c- n -Ip gauche NH17b — 2-me,3-aminopropene, c= c — c — n cis, c- c- n — Ip gauche NH18a — ethylenediamine, n — c- c — n anti, c-c-n-lpg+, g+ NH18b — ethylenediamine, n — c — c — n g +, c-c-n-lpg+,g+ NHl9a — N-methylaniline, N puckered NH19t — N-methylaniline, N planar NH20a — methylethylamine N-oxide, c- n — c — c anti NH20b — methylethylamine N-oxide, c — n — c — c gauche

## NH21 a — methylethylhydroxylamine, c — n — c- c anti

## NH21 b — methylethylhydroxylamine, c — n -c — c gauche NH22a — ethylamine N-oxide, o — n — c — c gauche NH22b — ethylamine N-oxide, o — n — c — c anti NH23a — ethylhydroxylamine, o — n — c — c gauche NH23b — ethylhydroxylamine, o — n — c — c anti Alcohols

## OH01 a — methanol OH02a — ethanol, c — c — o — h gauche OH02b — ethanol, c — c — o — h anti

## TABLE 1.

(continued) OH02j — ethanol, c — c — o — h = 0° OH02k — ethanol, c — c — o — h = 30° OH021 — ethanol, c — c — o — h = 60° OH02m- ethanol, c — c — o — h = 90° OH02n -ethanol, c — c — o — h = 120° OHO2o-ethanol,c-c-o-h = 150° OH02p — ethanol, c — c — o — h = 180° OH03a -n-propanol, c-c-c-o a, c-c-o-h g OH03b — n-propanol, c — c — c — o g -, OH03c — n-propanol, c — c — c — o g +, OH03d — n-propanol, c — c — c — o a, OH03e -n-propanol,c-c-c-og, c-c-0-h a OH04a — isopropanol, h — c — o — h gauche OH04b — isopropanol, h — c — o — h anti OH041 — isopropanol, h — c — o — h = 0° OH04k — isopropanol, h — c — o — h = 30° OH041 — isopropanol, h — c — o — h = 60° OH04m- isopropanol, h — c — o — h = 90° OH04n — isopropanol, h — c — o — h = 120° OH040 — isopropanol, h — c — o — h = 150° OH04p — isopropanol, h — c — o — h = 180° OH05a — t-butanol OH06a — cyclopentanol, equatorial Cs OH06b — cyclopentanol, axial Cs OH06c — cyclopentanol, equatorial C1 OH06d — cyclopentanol, axial C1 OH06j — cyclopentanol, Cs, ring planar OH07a — cyclohexanol, equatorial Cs OH07b — cyclohexanol, axial Cs OH07c — cyclohexanol, equatorial C1 OH07d — cyclohexanol, axial C1 OH08a — phenol OH09a -water OH1 Oa -vinyl alcohol, c= c — o — h trans OH1 Ob -vinyl alcohol, c= c — o — h skew OH1 l a — benzyl alcohol

## OH1 1 b — benzyl alcohol, h — o — c — c anti

## OH1 2a — propen-3-01, c= c — c — o skew,

## OH1 2b — propen-3-01, c= c- c- o cis,

## OH1 2c — propen-3-01, c= c — c — o skew,

## OH1 3a — 2-me-propen-3-01, c= c — c — o s, OH13b -2-me-propen-3-01, c= c-c-o c, OH14a — sec-butanol, ga I agb OH14b — sec-butanol, ga I ga

## OH1 4c — sec-butanol, ga I gg c-C-O-hg+ c-C-O-hg+ c-c-0-ha c-c-0-ha c-c-0-ha C-c-o-hg c-c-0-ha c-c-0-ha (Continues on next page)

## JOURNAL OF COMPUTATIONAL CHEMISTRY

501

## HALGREN

## TABLE I.

(continued)

## TABLE I.

(continued)

## OH1 4d — sec-butanol, ag I ag

## OH1 4e — sec-butanol, ag I ga

## OH1 4f — sec-butanol, ag I gg

## OH1 49 — sec-butanol, gg I ag

## OH1 4h — sec-butanol, gg I ga OH14i — sec-butaone, gg I gg

## OH1 4r — sec-butanol, cm I ag, approx ts OH14s -sec-butanol, cm Igg, approx ts

## OH1 4t — sec-butanol, ga I cm, approx ts

## OH1 4u — sec-butanol, ga I mp, approx ts

## OH1 4v — sec-butanol, ga I pc, approx ts OH14w-sec-butano1, mplag, approx ts

## OH1 4x — sec-butanol, mp I gg, approx ts

## OH1 4y — sec-butanol, pc I ag, approx ts

## OH1 42 — sec-butanol, pc 1 gg, approx ts OH15a — 1,2-ethanediol, h — o — c — c a, OH15b-1,2-ethanediol, h-0-c-cg-, OH15c — 1,2-ethanedioI_i h — o — c — c g -,

## OH1 5d — 1,2-ethanediol, h — o — c — c g+, Ethers

## OR01 a — methyl ethyl ether, c — o — c — c anti

## OR01 b — methyl ethyl ether, c — o — c — c gauche OR02a — methyl ethyl ether, c= c — o — c cis OR02b — methyl ethyl ether, c= c — o — c skew OR03a — diethyl ether, c — c — o — c anti,

## c -0 — c- c anti OR03b — diethyl ether, c — c — o — c anti, c-0-c-c gauche OR04a — methoxycyclohexane, equatorial Cs OR04b — methoxycyclohexane, axial C1 OR04c — methoxycyclohexane, equatorial C1 OR05a — oxetane, C2 OR05t — oxetane, planar OR06a — dimethyl ether OR07a — tetrahydrofuran, C2 OR07t — tetrahydrofuran, ring planar

## OR1 1 a — dioxolane, C2

## OR1 1 t — dioxolane, ring planar OR13a — methyl isopropyl ether,

## OR1 3b — methyl isopropyl ether, OR14a — methyl phenyl ether, c — o — c — c cis

## OR1 4j — methyl phenyl ether, c — o — c — c = 90° Alkanes RAOl a — methane RA02a -ethane, staggered RA02t — ethane, eclipsed 0-c-c-oa, c-c-0-h a 0-c-c-og+,c-c-0-ha o — c — c — o ~ +,c-C-o-hg+

## O-C-C-Og-,c-C-o-hg+ h — c — o — ch, gauche h-c-0-ch, anti ~ RA03a — propane RA04a — butane, c — c — c — c anti RA04b — butane, c — c — c — c gauche RA04t — butane, c — c — c — c = 0° RA05a — isobutane

## RAO6a — cyclobutane

## RAO6t — cyclobutane, ring planar RA07a — cyclopentane, half-chair C2 RA07t — cyclopentane, ring planar RA08a — cyclohexane, chair RA08b — cyclohexane, twist-boat C2

## RA1 Oa — methylcyclohexane, equatorial RAl Ob — methylcyclohexane, axial

## RA11 a — neopentane RA12a — 2,3-dimethylbutane, RAl2b — 2,3-dimethylbutane, h — c2 — c3 — h anti RA13a — cyclopropane RA14a — cyclooctane, crown D4d RA14b — cyclooctane, boat-chair Cs RA14c — cyclooctane, twist-boat-chair C2 RA14d — cyclooctane, S4 RA15a — methylcyclobutane, equatorial RA15b — methylcyclobutane, axial RA15j — methylcyclobutane, ring planar RA16a — cyclononane, [144] C2 RA16b — cyclononane, [3331 D3 RAl6c — cyclononane, [2251 C2 RA16d — cyclononane, [2341 C1 RA16e — cyclononane, [gal C1 Alkenes

## RE01 a — ethylene RE02a — propene RE03a — 1 -butene, c= c — c — c cis RE03b — 1 -butene, c= c — c — c skew RE04a — 1 -pentene, c — c — c — c anti RE04b — 1 -pentene, c — c — c — c gauche RE05a — 2-methyl-1 -butene, c= c — c — c skew RE05b — 2-methyl-1 -butene, c= c — c — c cis

## REO6a — isobutene RE07a — 1,Cpentadiene, c= c — c — c s +, s — RE07b — 1,4-pentadiene, c= c — c — c s -, s — RE08a — trans-2-butene RE08b — cis-2-butene

## REO9a — cyclobutene

## RE1 Oa — trans-2-pentene, c — c — c= c skew

## RE1 Ob — cis-2-pentene, c — c — c= c skew

## RE1 1 a — 1,4-~yclohexadiene Thiols, sulfides, and disulfides SROl a — hydrogen sulfide SR02a — ethanethiol, c — c — s — h gauche h — c2 — c3 — h gauche (used in validating MMFF94) (Continues on next page) 502

## VOL. 17, NOS. 5 & 6

## MERCK MOLECULAR FORCE FIELD. I

## TABLE 1.

(continued) SR02b — ethanethiol, c — c — s — h anti SR03a — dimethyl sulfide SR04a -ethyl methyl disulfide, c-c-s-s gauche SR04b — ethyl methyl disulfide, c — c — s — s anti SR05a — methyl hydrogen disulfide

## SRO6a — dimethyl disulfide SR07a — thiophenol SR07t — thiophenol, planar SR08a — methyl phenyl sulfide, c — c — s — c ca. 90°” SR08t — methyl phenyl sulfide, c — c — s — c = 0”

## SRO9a — 1 -propanethiol, s — c — c — c anti,

## SRO9b — 1 -propanethiol, s — c — c — c g -,

## SRO9c — 1 -propanethiol, s — c — c — c gauche,

## SR1 Oa — methyl hydrogen sulfide SRll a — 1,2-ethanedithioI_i all anti SRll b — 1,2-ethanedithiol, h — s — c — c anti, anti; s — c — c — s gauche SRll c — 1,2-ethanedithioI_i h — s — c — c anti, g +; s-c-c-sg- SR12a — methyl propyl sulfide, c — s — c — c gauche, s — c- c — c anti SR12b — methyl propyl sulfide, c — s — c — c g -, s-c-c-cg- SR12c — methyl propyl sulfide, c — s — c — c anti, s-c-c-c gauche SR12f — methyl propyl sulfide, c — s — c — c anti, s — c — c — c anti SR13a — thiomethanol h -s-c-c gauche h-S-C-Cg+ h — s-c-c anti aFor brevity, the conformational abbreviations a, g, t, c, and s are sometimes used for anti, gauche, trans, cis, and skew, respectively. These designations correspond approximately to torsion angles of 180”, 60”, 180”, 0”, and 120”. Where appropriate, the abbreviations “g+ ” and “g — ” or “s+ ” and “s — ” are used to indicate the relative signs of gauche or skew angles. bFor the OH14 conformers (sec-butanol), in the conforma- tional designations “wx / yz,” “w” indicates the conforma- tion of the c — c- c — o angle, “x” that of the c — c — c — c angle, “y” that of the h-0-c-ch, angle, and “2” that of the h-0-c-ch, angle; the designations “m” and “p” conote angles of approximately — 120” and 120”, re- spectively. cations, amine N-oxides, hydroxylamines, hydrox- amic acids, amidines, guanidines, amidinium cations, guanidinium cations, imadazolium cations, aromatic hydrocarbons, and heteroaromatic com- pounds (cf. Table 11). As can be inferred from the listing in Table I, the structural coverage is quite broad for some of these chemical families but is limited for others. Many of the bifunctional com- pounds included in the parameterization are un- saturated analogs of families listed above, that is, conjugated alkenes and aromatic hydrocar- bons (e.g., styrenes); a, punsaturated variants of amides, imines, aldehydes, ketones, carboxylic acids, esters, and carboxylate anions; vinylic ethers, alcohols, amines, and esters; and allylic aldehydes, ketones, amines, and alcohols. Other bifunctional compounds include: P-ketoacids; P-hydroxyesters; dicarboxylic acids; 1,2-diols, 1,2-diamines, and 1,2-dithiols; and nonconjugated dienes. Å limited selection of alkanes, amines, ketones, halides, es- ters and ethers containing four- or five-membered rings has also been studied. Compounds contain- ing SO, and oxyphosphorus groups have been treated as a part of the extension of the core parameterization described in part V.27 The number of chemical families treated in the core parameterization of MMFF94 is therefore large-certainly over 20-and many, though by no means all, combinations of functional groups of interest to medicinal and chemical industry chemists have been treated. Nevertheless, an in- crease of severalfold in the number of core MMFF94 parameters would probably be needed to allow the core force field to handle virtually any organic compound of pharamaceutical interest. To make MMFF94 as useful as possible, we have extended the core force field: (i) by parameterizing MMFF against a large set of experimental struc- tures extracted from the Cambridge Structural Database and against additional computational data; and (ii) by implementing a well-defined set of default-parameter assignments and carefully calibrated empirical rules for parameters not de- fined by either the structural data or the additional computational data.27 Some indication of the resul- tant range of chemical structures covered by MMFF94 can be gleaned from an examination of Table 111, which characterizes the current MMFF atom types. We hope to broaden the core, computationally derived, parameterization of MMFF in future work. Even now, however, the current set of core param- eters is, we believe, significandy broader than is provided in other specifically parameterized force field^.^,^,'^ In addition, we believe that the breadth and quality of the extended parameterization com- pares favorably, for organic compounds, with that provided by other force fields that employ generic

## JOURNAL OF COMPUTATIONAL CHEMISTRY

503

## HALGREN

## TABLE II.

Classes of Compounds Included in MMFF94's Core Parameterization. Alkanes, alkenes, aromatic hydrocarbons, conjugated alkenes and aromatics, nonconjugated alkenes. Five- and six-membered heteroaromatics. Alcohols, phenols, ethers, 1,2-diols, vinylic. Amines, imines, vinylic amines, allylic amines alcohols and ethers, allylic alcohols a,p-unsaturated imines, amidines, guanidines, 1,2-diamines. Hydroxyl amines, amine N-oxides. Amides, dipeptides, ureas, imides,. Aldehydes, ketones, a,p-unsaturated aldehydes. Ketals, acetals, hemiketals, hemiacetals. Carboxylic acids and esters, vinylic esters, a, p-unsaturated amides, hydroxamic acids and ketones, allylic aldehydes and ketones a,@-unsaturated acids and esters, dicarboxylic acids. p-ketoacids, p-hydroxyesters. Thiols, sulfides, disulfides, 1,2-dithiols. Halides (chlorides and fluorides). Amine, imine, amidine, guanidine, pyridine,. Carboxylate anions, a,@-unsaturated. Various saturated and unsaturated four- and and imidazole cations carboxylate anions five-membered ring systems parameters or empirical rules to achieve broad nominal coverage.", Computational Data Used in Parameterizing MMFF94 The computational data employed in parame- terizing the core force field fall into five main categories5': 1. Calculations at the HF/6-31G* level4' for — 500 HF/6-31G*-optimized geometries. These calculations, carried out using Gauss- ian 88,51 Gaussian 90,52 or Gaussian 92,53 covered nearly all of the molecular structures and conformations listed in Table I and also included — 70 hydrogen-bonded dimers used in the parameterization of intermolecular in- teractions." 2. Calculations with full geometry optimization at the MP2/6-31G* level for -360 equilib- r_{ium} conformers. This level of theory has been shown to give geometries for standard organic functional groups that rival experi- ment in ac~uracy.'~,~~ Theoretical geometries are particularly suitable for use in force-field parameterization because they do not entail the assumptions and artificial restrictions sometimes made in deriving experimental geometries," do not require sometimes ill- defined corrections for effects of thermal mo- tion," and are unlikely to manifest the large errors to which experimental determinations occasionally are subject." 3. Calculations for — 380 MP2/6-31G*-opti- mized geometries carried out at the MP2/TZP level using triple-zeta plus polar- ization basis sets, and at the MP2 and MP4SDQ levels using modified 6-31G* basis sets. As described in part IV,26 the MP2/TZP calculations and the MP3 and MP4SDQ cor- rections obtained using the modified 6-31G* basis set were combined to yield composite energies that we refer to as "MP4SDQ/TZP" energies, where the quotation marks indicate an approximation to full MP4SDQ/TZP results. 4. Single-point MP2/TZP calculations carried out at — 1450 torsionally incremented ge- ometries, derived from MP2/6-31G* geome- tries and partially optimized using refined but not yet final MMFF94 parameters. 5. Very large basis set calculations on inter- molecular interactions in nonpolar systems obtained using highly correlated wavefunc- tions. The use of these data in the derivation of MMFF94 is described in what follows. Methodology Used in Parameterizing

## MMFF94

The MMFF94 energy expression presented in eq. (1) contains seven terms. For the five valence- coordinate terms, MMFF94 employs quadratic force constants for bond stretching, angle bending, stretch-bend interaction, and out-of-plane bend- 504

## VOL. 17, NOS. 5 & 6

## MERCK MOLECULAR FORCE FIELD. I

ing; reference values for bond stretching and angle
bending; and one or more of the V_{ij} V_{ij} and V, constants for torsion interactions. Grouping the V, terms together, MMFF94 therefore utilizes seven classes of valence-coordinate parameters. MMFF94 also employs bond-charge increments oKl in eq. (14) and atomic polarizabilities al in eqs. (9) and (12) to generate quantities used in evaluating non- bonded interactions. In all, then, MMFF94 employs nine classes of force-field parameters. This section outlines the approaches used to derive each such class of parameters, and specifies how the individ- ual approaches were combined to yield a mutually consistent set of parameters. Full details are given elsewhere.’~ -27

## NONBONDED VAN DER WALLS AND

## ELECTROSTATIC PARAMETERS

As noted earlier, representative values for the atomic polarizabilities a; and for the derived mini- mum-energy separations RT, for nonhydrogen atoms have previously been de~cribed.,~ A prelim- inary listing of the associated MMFF atom types has also been given;39 the current set is specified in Table 111. For aliphatic hydrogens, the (Y and R* parameters were determined by fitting to high- quality ab initio data on intermolecular interactions for the methanes5 and hydrogen dimers. The vdW parameters for the polar hydrogens in water were determined by requiring that the water dimer be described in geometric and energetic terms similar to those found in successfully employed water models. This vdW representation was then trans- ferred to other types of polar hydrogen atoms. Initial values for the bond-charge-increment pa- rameters oK l in eq. (14) were obtained by employ- ing the Biosym Consortium program PROBE35 to fit the x, y, and z components of the molecular dipole moments to quantities computed at the HF/6-31G* level and scaled, as discussed in part II?* by a factor of 1.10. To make the fit well determined, partial atomic charges for polar hy- drogens and other key terminal atoms (and hence the associated wKl parameters) were fixed at val- ues representative of ESP-fit (electrostatic potential fit)% or Mulliken charges obtained from the ab initio calculations. Crucially, final vdW and elec- trostatic parameters for atom types involved in hydrogen-bonding interactions were obtained by adjusting the initial values to fit appropriately scaled41 interaction energies and hydrogen-bond- ing geometries computed at the HF/6-31G*

## TABLE 111.


$$
MMFF94 Symbolic and Numeric Atom Types Atom type Symbolic Numeric Definition [coordination numberla {formal chargeIb CR c=c
$$


## CSP2

c=o

## C=N

## CGD

## C=OR

## C=ON

coo

## COON

coo0 c=os c=s

## C=SN

cs02 cs=o css

## C=P

## CSP

=C= HC

## HSI

## -0 — OR oc=o oc=c

## OC=N

oc=s ON02

## ON=O


$$
OS03 os02 1 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 5 5 6 6 6 6 6 6 6 6 6 6 Alkyl carbon [41 Vinylic carbon I31. Generic sp’ carbon [31 Generic carbonyl carbon [31 Imine-type carbon [31 Guanidine carbon [31 Ketone or aldehyde carbonyl carbon [31 Amide carbonyl carbon 131 Carboxylic acid or ester carbonyl carbon [31 Carbamate carbonyl carbon [31 Carbonic acid or ester carbonyl carbon [31 Thioester carbonyl carbon, double bonded to 0 131 Thioester carbon, double bonded to S [31 Thioamide carbon, double bonded to S 131 Carbon in > C = SO, 131 Sulfinyl carbon in >C = S = 0 131 Thiocarboxylic acid or ester carbon [31 Carbon doubly bonded to P 131 Acetylenic carbon [2l Allenic carbon 121 Hydrogen attached to Hydrogen attached to Generic divalent oxygen [2] Ether oxygen [2l Carboxylic acid or ester oxygen 121 Enolic or phenolic oxygen [21 Oxygen in 4 — C = N — moiety [21 Divalent oxygen in thioacid or ester [21 Divalent nitrate “ether” oxygen [21 Divalent nitrate “ether” oxygen 121 Divalent oxygen in sulfate group [21 Divalent oxygen in sulfite group 121 carbon 111 silicon [I] (Continues on next page)
$$


## JOURNAL OF COMPUTATIONAL CHEMISTRY

505

## HALGREN

## TABLE 111.

(continued) Atom type Symbolic Numeric Definition [coordination numberla {formal chargeIb

## TABLE 111.

(continued) Atom type Symbolic Numeric Definition [coordination numberla {formal chargeIb oso os=o -0s OP03 OP02

## OPO

## -0P

o=c

## O=CN

## O=CR

o=co

## O=N

o=s o=s= NR

## N=C

## N=N

## NC=O

## NC=S

## NN=C

## NN=N

F CI Br I S s=c s=o

## > S=N

so2 S02N

## SO3

## SO4

6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 9 9 10 10 10 10 11 12 13 14 15 16 17 17 18 18 18 18 One of two divalent oxygens Divalent oxygen in Other divalent oxygen attached to sulfur 121 Divalent oxygen in phosphate group [2l Divalent oxygen in phosphite group 121 Divalent oxygen, one of two oxygens attached to P 121 Other divalent oxygen attached to phosphorus [21 Generic carbonyl oxygen [ l I Carbonyl oxygen in amides [ I I Carbonyl oxygen in aldehydes and ketones [I I Carbonyl oxygen in acids and esters [ I I Nitroso oxygen [ l ] Doubly bonded sulfoxide oxygen 111 O=S on sulfur doubly bonded to, e.g., C [ l I Amine nitrogen [31 lmine nitrogen [21 Azo-group nitrogen [2] Amide nitrogen [31 Thioamide nitrogen 131 Nitrogen in N — N=C moiety Nitrogen in N — N = N moiety Fluorine [ l I Chlorine [ l ] Bromine [ l ] Iodine [ I 1 Thiol, sulfide, or disulfide Sulfur doubly bonded to Sulfoxide sulfur [31 Tricoordinate sulfur doubly bonded to N [31 Sulfone sulfur [41 Sulfonamide sulfur [41 Sulfonate group sulfur [41 Sulfate group sulfur [41 attached to sulfur [2] R(R0)S = 0 [2] with deloc. Ip 131 with deloc. Ip 131 sulfur [21 carbon [ l I = so2

## SNO

SI

## CR4R

## HOR

HO

## CR3R

## HNR

## HPYL

## H3N

## HNOX

## HOCO

## HOP

## PO4

## PO3

## PO2

PO

## PTET

P

## HN=C

## HN=N

## HNCO

## HNCS

## HNCC

## HNCN

## HNNC

## HNNN

18 18 19 20 21 21 22 23 23 23 23 24 24 25 25 25 25 25 26 27 27 28 28 28 28 28 28 Sulfone sulfur, doubly bonded Sulfur in nitrogen analog of Silicon 141 Aliphatic carbon in 4-membered ring 141 Hydroxyl hydrogen in alcohols [I1 Generic hydroxyl hydrogen 11 I Aliphatic carbon in 3-membered ring 141 Generic hydrogen on sp3 nitrogen, e.g., in amines [I] Hydrogen on nitrogen in pyrrole [ l I Hydrogen in ammonia [ I 1 Hydrogen on N in a N-oxide Hydroxyl hydrogen in carboxylic acids [ l ] Hydroxyl hydrogen in

## H — 0 — P moiety [ 11 Phosphate group phosphorus [41 Phosphorus with 3 attached oxygens [41 Phosphorus with 2 attached oxygens [41 Phosphine oxide phosphorus [41 General tetracoordinate phosphorus [41 Phosphorus in phosphines 131 Hydrogen on imine nitrogen [I1 Hydrogen on azo nitrogen [ l l Hydrogen on amide nitrogen [A] Hydrogen on thioamide nitrogen [I1 Hydrogen on enamine nitrogen [ I ] Hydrogen in H — ## N — ## C= N

moiety [ I I Hydrogen in H — ## N — ## N=C

moiety [ l 1 Hydrogen in H — ## N — ## N = N

moiety [1 I to carbon a sulfone (Continues on next page 506

## VOL. 17, NOS. 5 & 6

## MERCK MOLECULAR FORCE FIELD. I

## HNSO

## HNC%

## HSP2

## HOCC

## HOCN

## CR4E

## HOH

02CM

## ONX

## O=N

02N 02N0 03N 0-s 02s 03s 04s

## OSMS

OP 02P 03P 04P 04CI

## HOS

## NR +

OM

## OM2

## HNR+

28 28 28 29 29 30 31 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 33 34 35 35 36 Hydrogen on NSO, NS02, or NS03 nitrogen 111 Hydrogen on N triply bonded to c [ l l Generic hydrogen on sp nitrogen [I] Enolic or phenolic hydroxyl hydrogen 111 Hydroxyl hydrogen in

## HO-C=N

moiety [ l l Olefinic carbon in 4-membered ring 131 Hydroxyl hydrogen in water 111 Oxygen in carboxylate group [11-1 / 2) Oxygen in N-oxides [ l I Oxygen in nitroso group [ I 1 Oxygen in nitro group [ l 1 Nitro-group oxygen in Nitrate anion oxygen [ l I Single terminal 0 on tetra- coordinate sulfur [ l I One of 2 terminal 0’s on sulfur [I I {variable)c One of 3 terminal 0’s on sulfur [ l l {variable)c Terminal 0 in sulfate anion [ l l (-1 /21 Terminal oxygen in thiosul-
finate anion 11 1 { — 1 I21 Oxygen in phosphine oxide [ l l One of 2 terminal 0’s on P 11 I {variableIc One of 3 terminal 0’s on P 11 1 {variable)c One of 4 terminal 0’s on P [ l 1 {variableIc Oxygen in perchlorate anion Hydrogen on oxygen attached
Quaternary nitrogen [41 {I 1 Oxide oxygen on sp3 carbon [ l I — 11 Oxide oxygen on sp2 carbon Hydrogen on quaternary nitrate 11 1 1-1 /31 [ l l (-1 /41 to sulfur [l 1 111 (-11 nitrogen [ I 1

## TABLE 111.

(continued) Atom type Symbolic Numeric Definition [coordination numberla {formal chargeIb

## TABLE 111.

(continued) Atom type Symbolic Numeric Definition [coordination numberla {formal chargeIb

## HIM +

## HPD+

## HNN +

## HNC +

## HGD +

CB

## NPYD

## NPYL

## NC=C

## NC=N

## NC=N

## NC%C

C02M

## CS2M

## NSP

NS02 NS03

## NC%N

## STH I

## NO2

## NO3

## N=O

## NAZT

## NSO

O+

## HO+

0=+

## HO= +

=N= 36 36 36 36 36 37 38 39 40 40 40 40 41 41 42 43 43 43 44 45 45 46 47 48 49 50 51 52 53

## N +=C

54 Hydrogen on imidazolium Hydrogen on pyridinium Hydrogen on amidinium Hydrogen on protonated Hydrogen on guanidinium Δromatic carbon, e.g., in Δromatic nitrogen with (T Δromatic 5-ring nitrogen with Enamine or aniline nitrogen, Nitrogen in N — C=N with Nitrogen in N — C = P with Nitrogen attached to C — C Carbon in carboxylate Carbon in thiocarboxylate Triply bonded nitrogen [ l I Sulfonamide nitrogen 131 Sulfonamide nitrogen 131 Nitrogen attached to cyano group 131 Δromatic 5-ring sulfur with T lone pair [21 Nitrogen in nitro group 131 Nitrogen in nitrate group [31 Nitrogen in nitroso group 121 Terminal nitrogen in azido or Divalent nitrogen replacing Oxonium oxygen [31 (1) Hydrogen on oxonium Oxenium oxygen [21 (1 1 Hydrogen on oxenium Central nitrogen in C=N = N
lminium nitroaen [3] { 1 } nitrogen 111 nitrogen [I1 nitrogen [ l l imine nitrogen 11 1 nitrogen [I1 benzene 131 lone pair [21 m- lone pair [21 deloc. Ip 131 deloc. Ip [31 deloc. Ip [31 triple bond 131 anion [31 anion [31 diazo group [ l I monovalent 0 in SO, group I21 oxygen [I1 oxygen [I1 or N=N=N 121 (Continues on next page)

## JOURNAL OF COMPUTATIONAL CHEMISTRY

507

## HALGREN

## TABLE 111.

(continued) Atom type Symbolic Numeric Definition [coordination numberla {formal chargeIb

## N +=N

## NCN +

## NGD+

## CGD+

## CNN +

## NPD +

## OFUR

## C% — ## NR%

NM

## C5A

## C5B

## N5A

## N5B

N20X N30X

## NPOX

## OH2

HS

## HS=N

HP s-P SM

## SSMO

S02M

## SSOM

=s=o

## -P=C

54 55 56 57 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 71 71 72 72 72 73 73 74 75 Positively charged nitrogen Either nitrogen in Guanidinium nitrogen Guanidinium carbon [31 Carbon in +N=C-N: resonance structures [31 Δromatic nitrogen in pyridinium [31 (1 1 Δromatic 5-ring oxygen with T lone pair 121 lsonitrile carbon [ l l Isonitrile nitrogen 121 Anionic divalent nitrogen Δromatic 5-ring C, (Y to Δromatic 5-ring C, p to Δromatic 5-ring N, (Y to Δromatic 5-ring N, p to sp2-hybridized N-oxide sp 3-hybridized N-oxide Pyridinium N-oxide Oxygen in water 121 Hydrogen attached to Hydrogen attached to > S = Hydrogen attached to Terminal sulfur bonded to Anionic terminal sulfur Terminal sulfur in thiosulfinate Sulfur in anionic sulfinate Tricoordinate sulfur in anionic Sulfinyl sulfur, e.g., in Phosphorus doubly bonded doubly bonded to N [31 (1 1 N+=C-N: 131 (1 /2) 131 11 /3) 121 { — 1) N:, 0:, or S: [31 N:, 0:, or S: [31 N:, 0:, or S: [2l N:, 0:, or S: [2l nitrogen [31 nitrogen [4] nitrogen [31 sulfur [I] sulfur doubly bonded to N [ l l phosphorus [ l I P [ l l [ l l (-11 group 111 (-1 /2) group 131 thiosulfinate group [31 c=s=o to C [31

## TABLE 111.

(continued) Atom type Symbolic Numeric Definition [coordination numberla {formal chargeIb

## N5M

CL04 c5 N5

## CIM +

## NIM +

## N5A +

## N5B +

N5 +

## N5AX

## N5BX

N50X FE+2 FE+3 F-

## CL — ## BR — ## LI +

## NA +

K+ ZN+2

## CA+ 2

CU+l cu+2 MG+2 76 77 78 79 80 81 81 81 81 82 82 82 87 88 89 90 91 92 93 94 95 96 97 98 99 Nitrogen in 5-ring aromatic Perchlorate anion chlorine [41 General carbon in 5-membered General nitrogen in
anion [2] {variable)c heteroaromatic ring [3] 5-mem bered heteroaromatic ring [21 Δromatic carbon between N’s in imidazolium 131 Δromatic nitrogen in imidazolium [31 (1 /2) Positive nitrogen in 5-ring alpha position [31 (1 I Positive nitrogen in 5-ring alpha position [31 (1 I Positive nitrogen in other 5-ring position [31 (1 1 N-oxide nitrogen in 5-ring alpha position [31 N-oxide nitrogen in 5-ring beta position 131 N-oxide nitrogen in other 5-ring position [31 Dipositive iron cation [O] (2) Tripositive iron cation [OI (3)
Floride anion [Ol { -1 1 Chloride anion [O] -1 1
Bromide anion [O] { -1 1 Lithium cation [O] (1) Sodium cation [OI (1 1 Potassium cation [OI (1 1 Dipositive zinc cation [OI (21 Dipositive calcium cation 101 (21 Monopositive copper cation 101 (1) Dipositive copper cation [OI 12) Dipositive magnesium cation [OI (2) ~ ~~~ aNumber of attached atoms. blnitial full or fractional charge, from which final MMFF94 partial atomic charges are obtained by adding contributions arising from the relative polarity of bonds involving attached atoms. ‘The formal charge is determined by dividing the net ionic charge on the SO, or PO, group among the equivalent terminal oxygens. 508

## VOL. 17, NOS. 5 &6

## MERCK MOLECULAR FORCE FIELD. I

## GEOMETRIC PARAMETERS

The reference bond lengths, Y:,, and bond an- gles, IY~~,, that appear in eqs. (2) and (3) were determined as follows. Given a trial set of MMFF parameters, optimized MMFF geometries for the molecules used in the parameterization were ob- tained from, and then systematically compared to, the reference ab initio ge0metries.5~ For each dis- tinct type of bond or angle (as determined by the MMFF atom types and the “bond-type” or ”angle-type” indexE), the average signed devia- tion between the MMFF and the ab initio bond lengths or angles was then determined and was used to adjust the trial reference value. The itera- tive procedure was initiated by setting the trial reference values equal to the average of the actual bond lengths or angles observed in the ab initio structures. As discussed in part III,E this approach had to be modified slightly to determine reference angles in small-ring compounds. This procedure was applied both to the MP2/6-31G*-optimized structures and to a similar set of HF/6-31G*-opti- mized geometries (ca. 350° structures in each case); the reference bond lengths and angles derived from fitting to the HF/6-31G* geometries were used in the fits to the HF/6-31G* first and second derivatives described in the following subsection.

## QUADRATIC FORCE CONSTANTS

Force constants for bond stretching, angle bend- ing, stretch-bend interaction, and out-of-plane bending were determined by using the Biosym Consortium program PROBE35 to fit a slightly modified version of the MMFF94 energy expres- sion to the Cartesian first and second derivatives of the HF/6-31G* energy. The principal modifica- tion consisted in replacing MMFF94’s Buf-14-7 and buffered electrostatic terms by Lennard-Jones
10-6 and simple coulombic [ 6 = 0 in eq. (1311 terms. As earlier work has suggested that valence- coordinate force constants are not strongly affected even by the neglect of nonbonded  interaction^:^ the substitution of comparable terms seems un- likely to have had an appreciable effect on the derived force constants. In these fits, only the quadratic force constants were optimized; parame- ters of all other classes were held constant. Finally, the HF/6-31G*-derived quadratic force constants were modified for use in MMFF94 by applying scaling factors chosen to optimize the fit of MMFF to experimental vibrational frequencies. Further details are given in part III.%

## TORSION PARAMETERS

The V_{ij} V_{ij} and V, parameters in eq. (7) were derived from fits to conformational energies using TORFIT?* These fits used ”penalty function” re- straints in connection with a ”build-up” protocol in which all but certain twofold parameters ini- tially were given zero values. The ab initio refer- ence data consisted of relative conformational en- ergies, nearly all of which were determined either from the composite ”MP4SDQ/TZP” calculations carried out at MP2/6-31G*-optimized geometries for — 380 conformers (Set A) or from single-point MP2/TZP calculations carried out at — 1450 tor- sionally incremented geometries derived from MP2/6-31G*-optimized geometries (Set B). Bench- mark calculations using still higher levels of the- ory and comparisons to experiment showed these to be the best tractable levels currently available to us.26 Set Å afforded 249 comparisons of “MP4SDQ/TZP” energies for optimized equilib-
r_{ium} or torsionally constrained conformers. Set B in turn yielded 1192 energy comparisons, each of which relates the MP2/TZP energy of a structure derived from a MP2/6-31G*-optimized equilib-
r_{ium} conformer to that of a ”torsion profile” struc- ture obtained by rotating one torsion bond by a specified extent (e.g., k 30°, k 60°,... 1. The inclu- sion of these comparisons assured that MMFF94 has a reasonable understanding of torsional pro- files and barriers. Full details are given in part Iv.26 We view the determination of torsion parame- ters as a particularly strong component in the development of core MMFF94. No other force field, to our knowledge, has employed so broad a range of comparably accurate data on conformational energies in its derivation.

## DETERMINATION OF MUTUALLY CONSISTENT

## MMFF PARAMETERS

Most force fields, have been derived using a “functional group” approach in which, for in- stance, “hydrocarbon parameters” are determined by fitting to data on alkanes and are then frozen. When alcohols and ethers, for example, are fit, only the parameters that arise from the newly introduced oxygen and polar hydrogen atom types need to be determined. This approach greatly sim-

## JOURNAL OF COMPUTATIONAL CHEMISTRY

509

## HALGREN

plifies the derivation of the force field but fails to allow for the possibility that correlations between parameters may veld values that fit the limited original data (e.g., on hydrocarbons) well but are poorly defined and/or are not optimal for describ- ing subsequent data (e.g., for hydrocarbon frag- ments in alcohols and ethers, etc.). A better strategy would be to determine all the force-field parameters simultaneously from the full set of experimental and/or computational data. Such an approach would ensure that any short- comings in the performance of the force field would be attributable to its form, or to the quality of the data used, rather than to its means of parameteri- zation. This approach is computationally impracti- cal at the present time. Fortunately, however, many classes of force-field parameters depend only weakly on others. For example, quadratic force constants change modestly when small changes are made in molecular geometries, and reference bond lengths and angles are insensitive to values employed for torsion parameters. This weak de- pendence allowed us to fashion a composite strat- egy that provided a computationally tractable approximation to the ideal of simultaneous deter- mination of all parameters. We implemented this strategy by carrying out between three and four interactions over the set of procedures described in the previous three subsections; for good measure, we also redetermined the nonbonded parameters for hydrogen-bonding interactions, as described earlier, before the final determination of the tor- sion parameters. This approach allowed each class of parameters to be determined in the context of successively refined values for parameters belong- ing to other classes. As a result, nearly all parame- ters derived by this set of procedures have been determined in a physically self-consistent fashi0n.5~ Performance of MM2X and MMFF94 This section summarizes MMFF94's ability to reproduce ab initio data used in its parameteriza- tion and also notes how well MM2X46 performs. Further details may be found in the accompanying studies.24-26

## MOLECULAR DIPOLE MOMENTS

For MMFF94, the partial atomic charges calcu- lated from the computationally derived bond charge increments reproduced the set of 423 HF/6-31G* molecular dipole moments, increased by 10% as described above, with the r_{ms} devia- tions6' shown:

## MMFF94

## MM2X

Dipole magnitude 0.39 D 0.64 D Dipole direction 5.5° 10.8° Also listed are the results obtained using the less widely parameterized MM2X force field for a somewhat smaller set of HF/6-31G* dipole mo- ments. For comparison, the r_{ms} value of the scaled HF/6-31G* dipole moments is 3.42 D. Thus, the average MMFF94 error is slightly larger than lo%, while that for MM2X is closer to 20%. MMFF94 is also considerably more accurate for dipole direc- tions. The present performance, and that for inter- molecular interaction energies and geometries in hydrogen-bonded dimers,24 appears quite reason- able for an approach that is simple enough to allow virtually automatically application to a wide range of organic and bio-organic systems. Never- theless, the treatment of electrostatic interactions is one area in which improvement particularly needs to be made in the future.

## EQUILIBRIUM BOND LENGTHS

The comparisons shown below involved a total of 4205 equilibrium bond lengths. They were ob- tained by using MMFF to optimize 358 MP2/ 6-31GY equilibrium conformers and by systema- tically comparing the ab initio- and MMFF- optimized geometries.25 For MM2X, 324 conform- ers covered by its parameterization were opti- mized, yielding a comparison of 3850 bond lengths. The results, cited below, are stated as r_{ms} devia- tions in angstroms from the optimized MP2/6-31G* bond lengths:

## MMFF94 MM2X

Equilibrium bond lengths 0.006 0.01 8 Clearly, the results for MMFF94 are excellent, those for MM2X respectable.

## EQUIIJBRIUM BOND ANGLES

A total of 7021 equilibrium bond angles for MMFF94 and 6462 bond angles for MM2X were
examined. The results are stated as r_{ms} deviations 51 0

## VOL. 17, NOS. 5 & 6

## MERCK MOLECULAR FORCE FIELD. I

in degrees from the optimized MP2/6-31G* bond angles:

## MMFF94 MMW

Equilibrium bond anQles 1.1 6° 1.70° Here, too, the MMFF94 results are quite good. MM2X also performs reasonably well.

## WILSON OUT-OF-PLANE PUCKERING ANGLES

For MMFF94, 237 conformers had out-of-plane centers involving a total of 1926 out-of-plane an- gles. For MM2X, 206 conformers had 1755 such angles. Many of the comparisons are of little inter- est, however, as all the methods find carbonyl and olefinic carbon to be essentially planar in un- strained compounds. In contrast, nonplanarity at nitrogen is found in the MP2/6-31G* structures for aliphatic amines, for most amides (which in this work include hydroxamic acids), and for such "unsaturated" amines as amidines, guanidines, vinylic amines and aromatic amines. For these
classes, the following MMFF94 and MM2X6' r_{ms} deviations were found:

## MP2

## MMFF94

r_{ms} angle r_{ms} dev. Δrnides (1 83 angles) 22.6° 9.38° Unsaturated arnines (33 angles) 43.8° 2.05° Saturated arnines (96 angles) 57.5° 0.91°

## MP2

## MM2X

r_{ms} angle r m s dev. Δrnides (1 53 angles) 18.5° 17.5° Unsaturated arnines Saturated arnines (33 angles) 43.8° 43.8° (84 angles) 56.9° 2.24° Shown for comparison are the r_{ms} values of the MP2 /6-31G*-op timized Wilson angles. Clearly, MMFF94 is far superior, though even it encounters some difficulty with amides, whose nitrogen center is notoriously easy to deform.62 As
we show in part III,= however, MMFF94 gives r_{ms} values for Wilson angles in primary, secondary, and tertiary amides that correctly reproduce the degree of puckering found in the MP2/6-31G* structures in an overall sense. MMFF94 also cor- rectly describes the nonplanar equilibrium geome- tries of unsaturated amines, whereas MM2X does not.

## TORSION ANGLES

Comparisons for a total of 7974 torsion angles for MMFF94 and 7409 for MM2X gave r_{ms} devia- tions from the MP2/6-31G*-optimized torsion an- gles as stated below:

## MMFF94

## MM2X

Torsion angles 5.83° 1 1.38°

$$
As discussed in part III,= MMFF94 sometimes underestimates and sometimes overestimates the degree of pyramidalization at nitrogen found in the MP2/6-31G* structures for amides. The asso- ciated errors in out-of-plane angles also affect the computed torsion angles and contribute signifi-
candy to the cited overall r_{ms} deviation. In addi- tion, in a number of cases involving methyl groups attached to sp2-hybridized centers, the relative en- ergies for the torsionally incremented (Set B) struc- tures on the MP2/TZP surface suggest that the MP2/6-31G* geometries (to which comparison is being made) are not equilibrium conformers on the higher level surface (from which the torsion parameters were largely derived). When question- able cases involving methyl rotations are excluded,
the r_{ms} deviation for MMFF94 falls to about 5°.26
$$


## CONFORMATIONAL ENERGIES AND TORSION

## PROFILES

As previous described, the torsion parameters were derived in fits to two sets of data on confor- mational energies. The results, stated as r_{ms} devia- tions in kilocalories per mole, were as follows:

## MMFF94

## MM2X

~~ Conforrnational energies Torsion-profile energies (Set A) 0.31 1.12 (Set 6) 0.50 1.57 For comparison, r_{ms} values for the relative ener- gies were 3.88 kcal/mol for the conformational energies and 4.37 kcal/mol for the torsion-profile energies for MMFF94, and 2.30 and 4.38 kcal/mol, respectively, for MM2XF3 Thus, MMFF94 accounts for about 90% of the variation in the ab ab initio

## JOURNAL OF COMPUTATIONAL CHEMISTRY

51 1

## HALGREN

relative energies in each case, whereas MM2X ac- counts for about 50%. We note that Sets Å and B each contained exten- sive conformational data on the glycine and ala- nine dipeptide analogs@ and on the full, methyl- capped glycine and alanine dipeptides. As we show in part IV,26 MMFF94 reproduces these data very well. All common protein sidechains are also cov- ered in its parameterization. No other published force field, to our knowledge, has been derived from a comparably extensive set of high-quality data on conformational comparisons pertinent to simulations on proteins.

## ADDITIONAL COMPARISONS

For the linear water dimer (optimized with the

## O-H... 0 angle restricted to 180°), MMFF94 gives a dimerization energy,Of -6.53 kcal/mol, an 0... 0 distance of 2.75 A, and an angle between the 0... 0 axis and the acceptor H-0-H plane of 27".24 The analoogous q ~ a n t i t i e s ~ ~ are — 6.50 kcal/ mol, 2.74 Å and 27°  for TIP3P water; -6.59 kcal/ mol, 2.75 Å and 26°  and for water; and -6.24 kcal/mol, 2.74 Å and 46°  and for TIP4P water. Thus, "MMFF water" behaves similarly in this static test. Work using liquid-phase simulations is currently underway to test and, if necessary, to reformulate or reparameterize MMFF94.66 Results for geometries and interaction energies for an extensive series of hydrogen-bonded dimers are presented in part II.', The comparisons show that MMFF94 accurately reflects the trends in in- teraction energies and geometries manifested in the ub initio calculations. The force field therefore appears to properly balance the strengths of wa- ter-water, water-solute, and solute-solute inter- actions. These comparisons suggest that MMFF94 can be used with confidence in computational studies of ligand-receptor binding. Also given in part 1124 are comparisons of vdW interaction ener- gies for the (CH,), and (H2), homodimers as a function of separation and orientation. These com- parisons show that MMFF94 accounts reasonably well for prototype nonpolar vdW interactions. Accuracy in Predicting Experiment This section summarizes MMFF94's ability to reproduce experimental data. Further details may be found in parts 111 and IV25,26; comparisons to experiment for the extended MMFF94 parameteri- zation are given in part V.27

## MOLECULAR GEOMETRIES

We have compared MMFF94 to experiment and to published MM3 geometries for a series of 30 organic molecules covering a variety of functional groups. For bond lengths, r_{ms} deviations reJative to experiment wereo found to be 0.014 Å for MMFF94 and 0.010 Å for MM3; for bond angles, the r_{ms} deviation was 1.2°  for each force field. Thus, MMFF94 is as successful as MM3 in predict- ing experimental bond angles, despite the fact that no experimental data on molecular geometries was used in deriving MMFF94. MM3 predicts experi- mental bond lengths more accurately, even in this test in which the experimental bond lengths were not strictly limited to the rg values MM3 seeks to emulate, but whether a difference of this magni- tude is of practical significance for molecular sim- ulations is unknown. Some of the difference in predicting experimental bond lengths arises from small, systematic deviations from experiment in the underlying MP2/6-31G* bond lengths. Part arises from the intrinsic difference between en- ergy-minimized (MP2/6-31G*) and thermally av- eraged (experimental) bond lengths; as a force field intended for use in molecular-dynamics sim- ulations, MMFF94 reproduces the energy-mini- mum bonds lengths obtained from the ab initio calculations, whereas MM3 incorporates thermal- averaging effects into its static model. For torsion angles, one particularly notable difference occurs for the "cisoid" conformation of 1,3-butadiene, for which MMFF94 predicts a nonplanar energy-mini- mized structure, wheeas MM3 gives a planar structure. Å more complete discussion, including comparisons to CHARMm9 and UFF" and more detailed comparisons to MM3, is given in part 111.~~

## VIBRATIONAL FREQUENCIES

To further characterize MMFF94's performance, vibrational frequencies were calculated for forma- mide, benzene, formic acid, formaldehyde, acetal- dehyde, methylamine, ammonia, methanol, water, methane, ethane, ethylene, hydrogen sulfide, guuche-ethanethiol, and dimethyl disulfide. When compared with published MM3 and experimen- tally determined freq~encies?~ r_{ms} deviations ver- sus experiment were found to be 61 cm-' for MMFF94,57 cm-' for MM3 for the slightly smaller 51 2

## VOL. 17, NOS. 5 & 6

## MERCK MOLECULAR FORCE FIELD. I

subset of molecules for which MM3 vibrational frequencies had been published, and 60 cm-' for MMFF94 for the same set of molecules used in assessing MM3. Thus, MMFF94 and MM3 perform comparably on an overall basis. For MM3, how- ever, we noted a number of instances in which its parameterization had employed experimental fre- quencies that differed significandy-by nearly 400 cm-I
for a B_{ij} stretching mode in benzene- from other published experimental val- ues that themselves had been shown to be compat- ible with theoretically calculated frequencies. Such instances illustrate one of the hazards of deriving a force field from experimental data: such data, at times, can contain large errors that then become a part of the derived force field.25

## COFORMATIONAL ENERGIES AND

## ROTATIONAL BARRIERS

Energy differences calculated using MMFF94 reproduce a diverse set of 37 experimentally de- termined gas-phase and solution conformational energies, enthalpies, and free energies (r_{ms} value 2.3 kcal/mol), with an r_{ms} deviation of 0.38 kcal/ mol, as opposed to 0.37 kcal/mol for both the supporting "MP4SDQ/TZP" calculations and for MM3.26 Moreover, MMFF94 reproduces 28 experi-
mentally determined rotational barriers (r_{ms} value 3.7 kcal/mol) with a r_{ms} deviation of 0.39 kcal/ mol. Importandy, these comparisons, and others discussed in part IV,26 demonstrate that fitting MMFF94 to high-quality theoretical data has si- multaneously conferred the ability to fit experi- ment. MMFF94 can be expected to perform equally well for the wide range of systems for which it has been parameterized but for which little or no ex- perimental data are available. Implementation of MMFF in OPTIMOL, CHARMm, and BatchMin In this section, we discuss some pertinent ele- ments related to the implementation of MMFF94 in OPTIMOL,46 the host molecular-mechanics plat- form for which MMFF94 and MM2X were devel- oped; the same elements also apply to the recent implementations of MMFF93 in CHARMm% and of MMFF94 in Bat~hMin.2~ In each of these implementations, the user (or the invoking modeling platform) simply repre- sents the subject molecule in language familiar to the organic chemist, that is, as a collection of atoms joined by single, double, or triple bonds, some atoms of which may have a nonzero formal charge. Δromatic systems may be supplied in any constituent Kekule form. The program then uses the supplied structural information to generate all additional information needed to carry out the calculation. It automatically determines the tor- sional "tree structure," perceives and classifies rings, defines symbolic atom types based on local connectivity, detects aromaticity, and creates ap- propriate lists of bond, angle, and torsional inter- actions. As previously des~ribed?~ the symbolic atom types (cf. Table 111) are then translated into the numeric values used to assign force-field pa- rameters to the force-field interaction terms. In establishing the relationship between param- eters and force-field interactions, the parameter files, which are kept in "canonical order" based on indices derived from the numerical atom types, are processed using a rapid binary search algo- rithm. If present, the fully qualified parameter cor- responding to the precise set of atom types (sup- plemented, in ambiguous situations, by defined bond, angle, stretch-bend, or torsion "interaction type^''*^^^^) is retrieved and used. For vdW_i bond stretching, stretch-bend, and bond-charge- increment parameters, no equivalences are recog- nized. For angle bending, out-of-plane bending, and torsion interactions, however, MMFF94 exe- cutes a staged "step-down" procedure in which increasingly generic values are sought whenever the "fully qualified" parameter is not found. This protocol is governed by the entries in Table IV, where the "Level 1°  atom types define the fully qualified parameters. Entries from Levels 2-5 are employed as needed in subsequent searches; those at Level 5, always "0°  except for atomic ions, serve as wild cards. Such wild card values are used only for peripheral atoms in an angle bending or torsional interaction or for noncentral atoms in an out-of-plane interaction. Level 4 generally corre- sponds to the atomic species, and Level 3 to atomic species plus hybridization. Currently, the first two levels employ identical numerical atom types. Unique values for Level 1 may be specified later if certain atom types need to be defined more specif- i ~ a l l y. ~ ~ The protocol used in the step-down proce- dure depends on the type of interaction (angle, torsion, out-of-plane).68 If no parameter is found, one of a series of carefully calibrated empirical rules is invoked (cf. part VZ7?. This staged- search/default-rule procedure allows applications to go forward when specific parameters are un-

## JOURNAL OF COMPUTATIONAL CHEMISTRY

51 3

## HALGREN

## TABLE IV.

Numerical Atom Type Equivalences Used in Assigning MMFF94 Parameters

## MMFF

Equivalence levelb symbola 1 2 3 4 5 CR c=c c=o

## CSP

HC OR o=c NR

## N=C

## NC=O

F CL BR I S s=c so s o 2 SI

## CR4R

## HOR

## CR3R

## HNR

## HOCO

## PO4

P

## HN=C

## HNCO

## HOCC

## CE4R

## HOH

02CM

## HOS

## NR +

OM

## HNR+

CB

## NPYD

## NPYL

## NC=C

C02M

## NSP

NS02 STHl

## NO2

## N=O

## NAZT

## NSO

O+

## HO+

0=+ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 1 21 22 23 24 25 26 28 28 29 2 31 7 21 8 6 36 2 9 10 10 3 42 10 16 10 9 42 9 6 21 7 1 1 1 1 5 6 6 8 8 8 11 12 13 14 15 15 15 15 19 1 5 1 5 5 25 25 5 5 5 1 31 6 5 8 6 5 1 8 8 8 1 8 8 15 8 8 8 8 6 5 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

## TABLE IV.

(continued) Equivalence levelb

## MMFF — sym bola 1 2 3 4 5

## HO= +

=N= N+=C

## NCN +

## NGD+

## CGD+

## NPD+

## OFUR

C%

## NR%

NM

## C5A

## C5B

## N5A

## N5B

N20X N30X

## N POX

## OH2

HS

## S2CM

S02M =s=o

## -P=C

## NM5

CL04 c 5 N5

## CIM +

## NIM +

## N5AX

FE+2 FE+3 F-

## CI — Br — ## LI +

K+

## NA +

ZN+2 CA+2

## CU+I

c u + 2 MG+2 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 87 88 89 90 91 92 93 94 95 96 97 98 99 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 87 88 89 90 91 92 93 94 95 96 97 98 99 21 42 9 10 10 2 10 6 4 42 10 2 2 9 9 9 8 9 70 5 16 18 17 26 9 12 2 9 2 10 9 87 88 89 90 91 92 93 94 95 96 97 98 99 5 8 8 8 8 1 8 6 1 8 8 1 1 8 8 8 8 8 70 5 15 15 15 25 8 12 1 8 1 8 8 87 88 89 90 91 92 93 94 95 96 97 98 99 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 70 0 0 0 0 0 0 0 0 0 0 0 0 87 88 89 90 91 92 93 94 95 96 97 98 99 ashown are representative MMFF94 symbolic atom types (cf. Table Ill). bThe Level 1 numerical atom types are the primary values. The usage of the equivalences reflected in Levels 2-5 is described in the “Implementation of MMFF94 in OPTIMOL, CHARMm, and BatchMin” section. 51 4

## VOL. 17, NOS. 5 & 6

## MERCK MOLECULAR FORCE FIELD. I

available, though inevitably with a loss in reliabil- ity. Concluding Discussion This and the accompanyig ~tudies’~-’~ intro- duce MMFF94, the initial published version of the Merck Molecular Force Field. As was noted in the Introduction, this version of MMFF is primarily intended for use in molecular-dynamics studies; a modified version intended for use in energy-mini- mization studies is under development?’ MMFF94‘s formulation and parameterization has a number of distinguishing features. One is that MMFF94 uses a unique functional form for describing van der Waals interactions and employs novel combination rules that embody a systematic correlation of vdW parameters with those that describe experimentally well-characterized interac- tions involving small molecules and rare-gas atoms.39 A second distinguishing feature is that MMFF94‘s core parameterization is primarily based on a large amount of computational data obtained from ab initio calculations-approximately 500 molecular structures optimized at the HF/6-31G* level, 475 structures optimized at the MP2/6-31G* level, 380 structures evaluated at the composite “MP4SDQ/TZP” levelz6 using MP2/6-3lG*-opti- mized geometries, and 1450 structures evaluated in single-point calculations at the MP2/TZP level. While ab initio data have been used in force devel- opment for at least two no other effort, to our knowledge, has used so much data of such high quality. A third distinguishing feature is that the core, computationally derived, portion of MMFF94 has been parameterized for an unusually wide variety of chemical systems. As a result, MMFF94 pro- vides well-defined parameters for more than 20 chemical families and treats many frequently oc- curring combinations of functional groups. The range of coverage for the extended parameteriza- tion is far larger stiiP7 The methodology used in parameterizing MMFF94 represents a fourth distinguishing fea- ture. Specifically, nearly all MMFF94 parameters have been determined in a mutually consistent fashion59 from the full set of available computa- tional data. Other force-field derivations have usu- ally employed a ”functional group” approach in which certain parameters are fit to a portion of the available data and are then frozen, While practical limitations of the ”functional group” approach have not yet been convincingly demonstrated, we prefer an approach that, by construction, yields mutually consistent values for the parameters. These attributes of its functional form and pa- rameterization combine to produce a force field that, by contemporary standards, performs very well. In particular, both computational data and experimental data are described well-the latter to a degree comparable to that achieved by MM3. These comparisons demonstrate that MMFF94’s parameterization against computational data has simultaneously conferred the ability to reproduce experiment. Consequently, MMFF94 can be ex- pected to perform equally well throughout the range of its parameterization from high-quality computational data, even for the many systems for which relevant experimental data is unavailable. This attribute constitutes a particularly strong ad- vantage of a computationally derived force field like MMFF94. Comparisons in functional form, performance, and/or manner of derivation for such other force fields as MM3, CFF93, OPLS, AMBER, CHARMm, UFF, and DREIDING are given in the accompanying st~dies.’~ -27 While the computational data used in its derivation necessarily relate to small molecules, it should be emphasized that MMFF94 has con- sciously been designed to be both a “small molecule” and a ”protein” force field. Among other factors, the excellent results obtained for conformational energies for dipeptide analogs and for dipeptidesZ6 and the uniform and balanced parameterization that has been pursued for non- bonded solvent-solvent, solvent-solute, and so- lute-solute  interaction^?^ when taken together with the reproduction of experimental data for small molecules with an accuracy comparable to that of MM3, suggest that MMFF94 should per- form well in both domains. Despite encouraging success, certain limitations are evident. One of particular importance arises from the fact that MMFF94 uses static atom- centered charges. As such, it neglects both higher order multipoles and electrostatic effects that arise from molecular polarizability. Because of these simplifications, MMFF94, like a number of other force fields, employs “enhanced” charge distribu- tions that emulate the effect of polarizability in amplifying electrostatic interactions for favorable contacts in a high-dielectric medium. Unfortu- nately, these enhanced charge distributions also amplify electrostatically unfavorable iteractions,

## JOURNAL OF COMPUTATIONAL CHEMISTRY

51 5

## HALG R EN

whereas proper account of polarizability would diminish them. They also improperly enhance electrostatic interactions in gas-phase or low-di- electric environments. Furthermore, they may not be optimal for describing intramolecular interac- tions, and may thereby limit the ability of the force field to account for differences in conformational energies. Indeed, compounds containing two or more strongly polar functional groups in close proximity have proven to be the most problematic in this respect, though good results have been obtained in most cases to date.z6 Other significant limitations include: the overly simplistic nature of the bond-charge-increment scheme used to assem- ble the partial atomic chargesz4; the lack of confor- mational dependence of the resultant chargesz4; and the omission of bond-torsion (and certain other) cross terms needed to account for significant geometrical changes that can occur when a torsion angle varies,25 an example being the elongation of ap amide partial C-N double bond by up to 0.1 Å when conjugation is broken. Å further signif- icant limitation is that no account is taken of metal-ligand interactions beyond that afforded by a relatively simplistic model that includes only electrostatic and van der Waals nonbonded interactions. What can be expected from future efforts at force-field development? First and foremost, better physical forms will need to be employed, particu- larly for electrostatic  interaction^.^" For example, even highly regarded water models such as SPC and TIP3P are known to describe certain configu- rations for the water dimer very poorly.71 In addi- tion, a broader selection of cross terms than are employed in MMFF94 will almost certainly be needed, and other enhancements can also be envi- si0ned.2~- 27 We expect that a computational ap- proach based almost solely on the use of ab initio data will become indispensable and that reliance on experimental data will diminish. The problem, ultimately, is one of information: too many force- field parameters, too little experimental data, and in many instances too nebulous a relationship be- tween the two. Fortunately, significant improve- ments in computer technology can be expected to make it increasingly feasible both to use more complex force fields in molecular simulations and to employ ever more rigorous computational mod- els to generate the data needed to derive them. But of course this approach will still yield a gas-phase force field, whereas most applications of interest to pharmaceutical and medicinal chemists take place in the condensed phase. This observation brings us back to an objective that underlies this work but has not yet been clearly stated: to define a forcefield that describes gas-phase molecular properties accurately and that behaves properly when the gas-phase system is embedded in the condensed phase. This objective can- not fully be met in a force field that treats electro- static interactions as simplistically as does the pre- sent version of MMFF. Ultimately, however, it will be met, because ”only” physics is involved, and because that physics is becoming increasingly well understood?’, Acknowledgment The author thanks the many colleagues at the Merck Research Laboratories who have supported this work; have contributed to the development of OPTIMOL, in which MMFF94 was first imple- mented; have developed computational tools used in MMFF94’s derivation; have created and main- tained the computational environment at Merck; or have utilized MMFF94 and its predecessor force fields in their own research. Many of these contri- butions are explicitly cited in the accompanying articles.z4 — z7 Supplementary Material Appendix Å (definition and role of the 16 MMFF94 parameter filesI7’ and Appendix B (com- puter-readable ASCII file containing the MMFF94 parameter filesz3 1 are available in Supplementary Material. References 1. J. B. Hendrickson, J. Am. Chem. Soc., 83, 4537-4547 (1961). For an even earlier reference, see: F. H. Westheimer, in Steric Effect in Organic Chemistry, M. S. Newman, Ed., Wiley, New York, 1956, Chapter 12. 2. D. H. Wertz and N. L. Allinger, Tetrahedron, 30, 1579 (1974). 3. (a) N. L. Allinger, \. Am. Chem. Soc., 89, 8127 (1977); (b) U. Burkert and N. L. Allinger, Molecular Mechanics; American Chemical Society, Washington, DC, 1982; (c) N. L. Allinger and Y. Yuh, QCPE, 12, 395 (1980). 4. N. L. Allinger, Y. H. Yuh, and J-H. Lii, 1. Am. Chem. SOC., 111, 8551-8566 (1989). See also: N. L. Allinger and L. Yan, J. Am. Chem. Soc., 115, 11918-11925 (1993), and references therein. 5. W. D. Cornell, P. Cieplak, C. I. Bayly, I. R. Gould, K. M. Merz Jr., D. M. Ferguson, D. C. Spellmeyer, T. Fox, J. W. Caldwell, and P. A. Kollman, I. Am. Chem. Soc., 117, 5179-5197 (1995); S. J. Weiner, P. A. Kollman, D. T. Nguyen, and D. A. Case, J. Comput. Chem., 7, 230-252 (1986); S. J. 51 6

## VOL. 17, NOS. 5&6

## MERCK MOLECULAR FORCE FIELD. I

6. 7. 8. 9. 10. 11. 12. 13. 14. 15. 16. 17. 18. 19. 20. 21. Weiner, P. A. Kollman, D. T. Nguyen, D. A. Case, U. C. Singh, C. Ghio, G. Alagona, S. Profeta Jr., and P. Weiner, (a) W. L. Jorgensen and J. Tirado-Rives, J. Am. Chem. SOC., 110, 1657-1666 (1988), and references therein; (b) H. A. Carlson, T. B. Nguyen, M. Orozco, and W. J. Jorgensen, J. Comput. Chem., 14, 1240-1249 (19931, and references therein. B. R. Brooks, R. E. Bruccoleri, B. D. Olafson, D. J. States, S. Swaminathan, and M. Karplus, J. Comput. Chem., 4, S. Lifson, A. T. Hagler, and P. Dauber, J. Am. Chem. SOC., 101,5111-5121 (19791, and references therein. F. M. Momany and R. Rone, J. Comput. Chem., 13, 888-900 (1992). S. L. Mayo, B. D. Olafson, and W. A. Goddard 111, J. Phys. Chem., 94, 8897 (1990).
A. K. RappC_i C. J. Casewit, K. S. Colwell, W. A. Goddard 111, and W. M. Skiff, J. Am. Chem. SOC., 114,10024-10035 (1992), and references therein. A. Vedani and D. W. Huhta, J. Am. Chem. SOC., 112, V. S. Allured, C. M. Kelly, and C. R. Landis, J. Am. Chem. SOC., 113, 1 (1991). D. M. R_{oot}, C. R. Landis, and T. Cleveland, J. Am. Chem. J. R. Maple, U. Dinur, and A. T. Hagler, Proc. Natl. Acad. Sci. USA, 85, 5350-5354 (1988). A. T. Hagler, J. R. Maple, T. S. Thacher, G. B. Fitzgerald, and U. Dinur, In Computer Simulation of Biomolecular Systems, W. F. van Gunsteren and P. K. Weiner, Eds., ESCOM, Leiden, 1989, pp. 149-167. U. Dinur and A. T. Hagler, in Reviews in Computational Chemistry, K. B. Lipkowitz and D. B. Boyd, Eds., VCH Publishers, New York, 1991, Vol. 2, pp. 99-163. See also: J. Palca, Nature, 322, 586 (1986). (a) J. R. Maple, M.-J. Hwang, T. P. Stockfish, U. Dinur, M. Waldman, C. S. Ewig, and A. T. Hagler, J. Comput. Chem., 15, 161-182 (1994); M.-J. Hwang, T. P. Stockfish, and A. T. Hagler, J. Am. Chem. SOC., 116, 2515-2525 (1994). N. L. Allinger, K. Chen, and J.-H. Lii, J. Comput. Chem. (this issue). (a) W. J. Hehre, L. Radom, P. v. R Schleyer, and J. A. Pople, Ab lnitio Molecular Orbital Theoy, Wiley, New York, 1986, Chapter 6; (b) D. J. DeFrees, B. A. Levi, S. K. Pollack, W. J. Hehre, J. S. Binkley, and J. A. Pople, J. Am. Chem. SOC., 101, For a recent example of a significant experimental error which was corrected on the basis of information from ab initio calculations, including calculations at the MP2/6-31G* level used in the present work, see: 8. J. Smith and L. Radom, J. Am. Chem. SOC., 112, 7525-7528 (1990). See, for example: (a) U. Dinur and A. T. Hagler, I. Chem. Phys., 91, 2949-2958 (1989); (b) U. Dinur and A. T. Hagler, 1. Chem. Phys., 91,2959-2970 (1989); (c) U. Dinur, J. Comput. Chem., 12,91-105 (1991); (d) U. Dinur, J. Comput. Chem., 12, 469-486 (1991); (e) U. Dinur, J. Phys. Chem., 94, 5669-5671 (1990). See, for example: (a) G. Corongiu, M. Migliore, and E. Clementi, J. Chem. Phys., 90, 4629 (1989); (b) L. X. Dang, J. E. Rice, J. Caldwell, and P. A. Kollman, J. Am. Chem. SOC., 113, 2481-2486 (19911, and references therein; (c) M. Sprik, J. Am. Chem. SOC., 106, 765-784 (1984). 187-217 (1983). 4759-4767 (1990). SOC., 115, 4201-4209 (1993). 4085-4089 (1979). 22. 23. 24. 25. 26. 27. 28. 29. 30. 31. 32. 33. 34. 35. 36. 37. 1. Phys. Chern., 95,2283-2291 (1991), and references therein. (d) S.-8. Zhu, S. Yao, J.-B. Zhu, S. Singh, and G. W. Robin- son, J. Phys. Chem., 95,6211-6217 (1991); (e) s. W. Rick, s. J. Stuart, and B. J. Beme, J. Chem. Phys., 101, 6141-6156 (1994); (f) D. N. Bemardo, Y. Ding, K. Krough-Jespersen, and R. M. Levy, J. Phys. Chem., 98, 4180-4187 (1994); (g) For an early implementation of induced-dipole effects, see also L. Dosen-Micovic, D. Jeremic, and N. L. Allinger, J. Am. Chem. Soc., 105, 1716, 1723 (1983). S. Dasgupta and W. A. Goddard 111, J. Chem. Phys., 90, The MMFF94 parameters (Appendix B, Supplementary Ma- terial) are available in computer-readable form (see foot- note * on first page of this article). Part 11: T. A. Halgren, J. Comput. Chem. (this issue). Part 111: T. A. Halgren, J. Comput. Chem. (this issue). Part IV T. A. Halgren and R. B. Nachbar, J. Comput. Chem. (this issue). Part V: T. A. Halgren, J. Comput. Chem. (this issue). This collaboration involved Prof. Martin Karplus (Harvard University) and Dr. Ryszard Czerminski and others of Molecular Simulations, Inc. (San Diego, CAI. Currently, a version of CHARMm that supports the earlier and less widely parameterized MMFF93 force field (which lacks, e.g., the ability to recognize a number of the ionic species parameterized in ref. 27; see also refs. 25 and 26) is avail- able from MSI. However, while the local Merck code for CHARMm employs MMFF94, arrangements for including MMFF94 in the distributed MSI version have not yet been concluded. P. S. Shenkin and T. A. Halgren (work in progress). The MacroModel program suite and its BatchMin module, de- veloped in the laboratories of Professor Clark Still, are available from Columbia University (New York, NY). OPTIMOL has been developed and maintained by the au- thor, but is based in part on computer code adapted from a public domain version of MM2 or written by Drs. R. B. Nachbar, B. L. Bush, G. M. Smith, E. F. Fluder Jr., and J. D. Andose of the Merck Research Laboratories. Distribution of OPTIMOL by the Quantum Chemistry Pro- gram Exchange (Indiana University) would permit free usage of the program but would prohibit its commercializa- tion. T. A. Halgren and R. 8. Nachbar (work in progress). The Merck Index, 11th ed., S. Budavari, Ed., Merck & Co., Rahway, NJ, 1989. Fine Chemicals Direct0 y Handbook, Fraser Williams (Scien- tific Systems), London, 1983-1985. Connection tables dis- tributed by Molecular Design Ltd., Hayward, CA. PROBE is a computer program used to derive molecular- mechanics parameters in least-squares fits to data obtained from ab initio calculations. PROBE was created for the Biosym Consortium on Potential Energy Functions by Biosym Technologies, Inc. (now Molecular Simulations, Inc.); cf. refs 15 and 16a. The derivation of MMFF94 used a 1991 version of PROBE. See, for example, G. C. Lie and E. Clementi, Phys. Rev., 33A, 2679 (1986). H. J. C. Berendsen, J. P. M. Postma, W. F. van Gunsteren, and J. Hermans, In Intermolecular Forces, B. Pullman, Ed., Reidel, Dordrecht, Holland, 1981, pp. 331-342. 7207-7215 (1989).

## JOURNAL OF COMPUTATIONAL CHEMISTRY

51 7

## HALGREN

38. B. M. Pettit, in Computer Simulation of Biomolecular Systems, W. F. van Gunsteren and P. K. Weiner, Eds., ESCOM, Leiden, 1989, pp. 94-100. 39. T. A. Halgren, J. Am. Chem. Soc., 114, 7827-7843 (1992). 40. W. J. Hehre, L. Radom, P. v. R Schleyer, and J. A. Pople, Å b Znitio Molecular Orbital Theory, Wiley, New York, 1986, Chapter 4. The 6-31G* basis sets also known as 6-31G(d). 41. A. D. MacKerell Jr. and M. Karplus, J. Phys. Chem., 95, 10559-10560 (1991); A. D. MacKerell Jr., J. Wi6rkiewicz- Kuczera, and M. Karplus, J. Am. Chem. Soc. (in press). 42. M. J. Frisch. J. E. Del Bene, J. S. Binkley, and H. F. Schaeffer 111, J. Chert]. Pkys., 84, 2279-2289 (1986). 43. T. A. Halgren, J. Am. Chem. Soc., 112, 4710-4723 (1990). 44. M. Orozco and F. J. Luque, 1. Comput. Chem., 14, 881-894 (1993). 45. E. B. Wilson Jr., J. C. Decius, and P. C. Cross, Molecular Vibrations, Dover, New York, 1955, Chapter 4. 46. M. K. Holloway, J. M. Wai, T. A. Halgren, P. M. D. Fitzger- ald, J. P. Vacca, B. D. Dorsey, R. B. Levin, W. J. Thompson, L. J. Chen, S. J. deSolms, N. Gaffin, A. K. Ghosh, E. A. Giuliani, S. L. Graham, J. P. Guare, R. W. Hungate, T. A. Lyle, W. M. Sanders, T. J. Tucker, M. Wiggins, C. M. Wiscount, 0. W. Woltersdorf, S. D. Young, P. L. Darke, and J. A. Zugay, J. Med. Chenz., 38, 305-317 (1995). 47. Note that this truncation does not lead to a "cubic-bend" catastrophe because the range of the angle is limited to 180° and the cubic-bend constant is relatively small. 48. See, for example, the "vdW and hydrogen bonding" pa- rameters listed in Table I of J.-H. Lii and N. L. Allinger, J. Comput. Chem., 12, 186-199 (1991). 49. M. Waldman and A. T. Hagler, J. Comput. Ckem., 14, 50. Most of the ab initio calculations used in parameterizing MMFF were performed on the Merck Research Laboratories Cray YPMSi/4-128 supercomputer. 51. M. J. Frisch, M. Head-Gordon, H. B. Schlegel, K. Rag- havachari, J. s. Binkley, C. Gonzalez, D. J. Defrees, D. J. Fox, R. A. Whiteside, R. Seeger, C. F. Melius, J. Baker, R. L. Martin, L. R. Kahn, J. J. P. Stewart, E. M. Fluder, S. Topiol, and J. A. Pople, Gaussian 88, Gaussian, Inc., Pittsburgh, PA, 1988, as modified at Merck for improved 1/0 performance by E. M. Fluder. 52. M. J. Frisch, M. Head-Gordon, G. W. Trucks, J. 8. Foresman, H. B. Schlegel, K. Raghavachari, M. Robb, J. S. Binkley, C. Gonzalez, D. J. Defrees, D. J. Fox, R. A. Whiteside, R. Seeger, C. F. Melius, J. Baker, R. L. Martin, L. R. Kahn, J. J. P. Stewart, S. Topiol, and J. A. Pople, GAUSSIAN 90 (Revision J), Gaussian, Inc., Pittsburgh, PA, 1990. 53. M. L. Frisch, G. W. Trucks, M. Head-Gordon, P. M. W. Gill, M. W. Wong, J. B. Foresman, B. G. Johnson, H. B. Schlegel, M. A. Robb, E. S. Replogle, R. Gomperts, J. L. Andres, K. Raghavachari, J. S. Binkley, C. Gonzalez, R. L. Martin, D. J. Fox, D. J. Defrees, J. Baker, J. J. P. Stewart, and J. A. Pople, GAUSSIAN 92 (Revision C), Gaussian, Inc., Pitts- burgh, PA, 1992. 54. Calculations at the MP2/6-31G*//MP2/6-31G* level are used as the basis of structure determination in the high-level composite G1 and G2 methods developed by Pople and coworkers; cf. L. A. Curtiss, K. Raghavachari, G. W. Trucks, and J. A. Pople, I. Chem. Phys., 94, 7221-7230 (1991), and references therein. Comparisons to experiment for a num- 1077-1084 (1993). 55. 56. 57. 58. 59. 60 61. ber of the molecules employed in the derivation of MMFF94 are summarized in A. St.-Amant, W. D. Cornell, T. A. Halgren, and P. A. Kollman, J. Comput. Chem., 16,1483-1506 (1995). See, for example, G. Chalasinski, M. M. Szczesniak, P. Cieplak,, and S. Scheiner, J. Chem. Phys., 94, 2873-2883 (19911, and references therein. The ESP-fit calculations were carried out with a version of Gaussian 88 to which Dr. M. D. Miller (Merck Research Laboratories) had interfaced PDM88 (D. E. Williams, QCPE, Program No. 568, 1988). In these MMFF optimizations, weak penalty-function re- straints were applied to the torsion angles to insure that comparable MMFF and nb initio conformations were being compared (cf. ref. 26). TORFIT is a versatile program developed at the Merck Research Laboratories which derives torsional parameters via least-squares fits to relative conformational energies (cf. ref. 26). The procedure used is not strictly "mathematically" self- consistent, however, because formal couplings between pa- rameters belonging to different classes (e.g., between refer- ence values and force constants for angles at trigonal cen- ters) have not been addressed. Further iterations would probably cause a slow drift away from the parameter val- ues reported in this work. We view the parameters as being "physically" self-consistent, however, in the sense that such further iterations would not materially improve the fit to the computational data.
The cited r_{ms} deviations in dipole directions are weighted
r_{ms} deviations constructed to avoid overemphasizing large errors in directions for dipole moments of small magnitude (cf. ref. 24). We should note that MM2X actually uses the Allinger (MM2/MM3) definition for out-of-plane angles. To clarify the comparison to MMFF, however, we have used the Wilson definition in analyzing the MM2X-optimized ge- ometries. The Allinger angles typically are about three times smaller in magnitude. 62. See, for example: M. W. Wong and K. B. Wiberg, J. Phys. Chem., 96, 668-671 (1992). 63. The first r_{ms} value is much lower for MM2X because some high-energy structures in the transition state region for

## C-N

amide bond rotation in N-methylformamine were poorly treated by MM2X and had to be removed from the test set (cf. ref. 26). 64. For structures, see T. Head-Gordon, M. Head-Gordon, M. J. Frisch, C. L. Brooks 111, and J. A. Pople, J. Am. Chem. Soc., 65. W. L. Jorgensen, J. Chandrasekhar and J. D. Madura, J. Chem. Phys., 79, 926-935 (1983). 66. (a) B. L. Bush, J. L. Banks, R. Czerminski, and T. A. Halgren (work in progress), using a local version of CHARh4m into which MMFF94 has been integrated; (b) T. A. Halgren, S.4. So, and M. Karplus (work in progress). 67. As an example, we might at some point wish to define different bond-charge increments for C=O groups in amides, esters, ketones, etc., for which differing symbolic atom types but common numeric atom types currently are assigned. The equivalence procedure provides a convenient way to do so without requiring that atom types and param- eters describing common bond, angle, torsion, and other interactions simultaneously be modified. 113, 5989-5997 (1991). 51 8

## VOL. 17, NOS. 5 & 6

## MERCK MOLECULAR FORCE FIELD. I

68. For bending of the i-j-k angle, a five-stage process based in the level combinations 1-1-1,2-2-2,3-2-3,4-2-4, and 5-2-5 is used. For i-j-k-2 torsion interactions, a five-stage process based on level combinations 1-1-1-1, 2-2-2-2, 3-2-2-5, 5-2-2-3, and 5-2-2-5 is used, where stages 3 and 4 correspond to "half-default" or "half-wild-card" en- tries. For out-of-plane bending ijk; I, where j is the central atom [cf. eq. (511, the five-stage protocol 1-1-1; 1, 2-2-2; 2, 3-2-3;3, 4-2-4;4, 5-2-5;5 is used. The final stage pro- vides wild-card defaults for all except the central atom. 69. N. L. Allinger and M. J. Hickey, Tetrahedron, 28, 2157-2161 (1972). Although not specifically referenced, the results of this work were used in developing a force field for carbonyl compounds: N. L. Allinger, M. T. Trible, and M. A. Miller, Tetrahedron, 28, 1173-1190 (1972) (N. L. Allinger, private communication). 70. W. F. van Gunsteren, In Computer Simulation of Biornolecular Systems, W. F. van Gunsteren and P. K. Weiner, Eds., ESCOM, Leiden, 1989, pp. 27-59. 71. M. Mezei and J. J. Dannenberg, 1. Phys. Chem., 92,5860-5861 (1988). 72. Appendix Å is available in Supplementary Material (see footnote * on the first page of this article).

## JOURNAL OF COMPUTATIONAL CHEMISTRY

51 9
