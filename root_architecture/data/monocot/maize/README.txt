This README.txt file was written on 20190225 by Benjamin Delory

-------------------
GENERAL INFORMATION
-------------------

1. Plant and root level data used for the parameterization and benchmarking of root architecture models simulating maize root systems.  

2. Author Information

        Name: Andrea Schnepf (principal investigator)
        Institution: IBG-3, Jülich Forschungszentrum GmbH
        Email: a.schnepf@fz-juelich.de

        Name: Guillaume Lobet (principal investigator)
        Institution: IBG-3, Jülich Forschungszentrum GmbH
        Email: g.lobet@fz-juelich.de

        Name: Benjamin M. Delory (postdoctoral research associate)
        Institution: Institute of Ecology, Leuphana University
        Address: Universitaetsallee 1, 21335 Lueneburg, Germany
        Email: benjamin.delory@leuphana.de

3. Source images: https://github.com/RSA-benchmarks/collaborative-comparison/tree/master/root_architecture/img/monocot/maize

4. Plant-level data were generated with the architect function of the R package archiDART (https://github.com/archidart/archidart) on 20190225

5. Root-level data were generated with the root function of the R package archiDART (https://github.com/archidart/archidart) on 20190225

-----------------------------------------------
DATA-SPECIFIC INFORMATION FOR: maize-plant.txt
-----------------------------------------------

1. Number of variables: 42

2. Number of cases/rows: 24

3. Column separator: ","

4. Variable List

    4.1. Name: FileName
         Description: the root image file name aggregated with the plant identification number (Ex: "PL01_DAS01_1" means plant 1 in the image file PL01_DAS01)
         Type: character

    4.2. Name: Time
         Description: the root system age
         Unit: days
         Type: integer

    4.3. Name: TRL
         Description: the total root system length
         Unit: cm
         Type: numeric

    4.4. Name: GRTR
         Description: the root system growth rate
         Unit: cm/day
         Type: numeric

    4.5. Name: L1R
         Description: the total first-order root length (primary root length + basal root length)
         Unit: cm
         Type: numeric

    4.6. Name: GR1R
         Description: the first-order root growth rate
         Unit: cm/day
         Type: numeric

    4.7. Name: TN1R
         Description: the number of first-order roots
         Type: integer

    4.8. Name: TNLR
         Description: the number of lateral roots
         Type: integer

    4.9. Name: TLRL
         Description: the total lateral root length
         Unit: cm
         Type: numeric

    4.10. Name: N2LR
          Description: the number of second-order roots
          Type: integer

    4.11. Name: L2LR
          Description: the total second-order root length
          Unit: cm
          Type: numeric

    4.12. Name: ML2LR
          Description: the average second-order root length
          Unit: cm
          Type: numeric

    4.13. Name: GR2L
          Description: the second-order root growth rate
          Unit: cm/day
          Type: numeric

    4.14. Name: D2LR
          Description: the lateral root density on the primary root
          Unit: root(s)/cm
          Type: numeric

    4.15. Name: Height
          Description: the root system height
          Unit: cm
          Type: numeric

    4.16. Name: Width
          Description: the root system width
          Unit: cm
          Type: numeric

    4.17. Name: Convexhull
          Description: the area of the convex hull
          Unit: cm^2
          Type: numeric

    4.18. Name: Magnitude
          Description: the magnitude of the root system (Fitter index)
          Type: integer

    4.19. Name: Altitude
          Description: the altitude of the root system (Fitter index)
          Type: integer

    4.20. Name: ExtPathLength
          Description: the total exterior path length of the root system (Fitter index)
          Type: integer

    4.21. Name: MD1
          Description: the average diameter of first-order roots
          Unit: cm
          Type: numeric

    4.22. Name: MD2
          Description: the average diameter of second-order roots
          Unit: cm
          Type: numeric

    4.23. Name: MDLR
          Description: the average lateral root diameter
          Unit: cm
          Type: numeric

    4.24. Name: S1
          Description: the total surface area of first-order roots
          Unit: cm^2
          Type: numeric

    4.25. Name: S2
          Description: the total surface area of second-order roots
          Unit: cm^2
          Type: numeric

    4.26. Name: Stot
          Description: the total root surface area
          Unit: cm^2
          Type: numeric

    4.27. Name: V1
          Description: the total volume of first-order roots
          Unit: cm^3
          Type: numeric

    4.28. Name: V2
          Description: the total volume of second-order roots
          Unit: cm^3
          Type: numeric

    4.29. Name: Vtot
          Description: the root system volume
          Unit: cm^3
          Type: numeric

    4.30. Name: LPR
          Description: the primary root length
          Unit: cm
          Type: numeric

    4.31. Name: LBR
          Description: the total basal root length
          Unit: cm
          Type: numeric

    4.32. Name: GRPR
          Description: the growth rate of the primary root
          Unit: cm/day
          Type: numeric

    4.33. Name: GRBR
          Description: the growth rate of the basal roots
          Unit: cm/day
          Type: numeric

    4.34. Name: TNBR
          Description: the number of basal roots
          Type: integer

    4.35. Name: TNPR
          Description: the primary root number (should be equal to 1)
          Type: integer

    4.36. Name: MLBR
          Description: the average length of basal roots
          Unit: cm
          Type: numeric

    4.37. Name: MDPR
          Description: the average diameter of the primary root
          Unit: cm
          Type: numeric

    4.38. Name: MDBR
          Description: the average diameter of the basal roots
          Unit: cm^3
          Type: numeric

    4.39. Name: SPR
          Description: the total surface area of the primary root
          Unit: cm^2
          Type: numeric

    4.40. Name: SBR
          Description: the total surface area of the basal roots
          Unit: cm^2
          Type: numeric

    4.41. Name: VPR
          Description: the volume of the primary roots
          Unit: cm^3
          Type: numeric

    4.42. Name: VBR
          Description: the total volume of basal roots
          Unit: cm^3
          Type: numeric

5. Missing data codes:
        NA        Not available

----------------------------------------------
DATA-SPECIFIC INFORMATION FOR: maize-root.txt
----------------------------------------------

1. Number of variables: 16

2. Number of cases/rows: 1914

3. Column separator: ","

4. Variable List

    4.1. Name: file
         Description: the root image file name
         Type: character

    4.2. Name: plant
         Description: the plant identification number
         Type: integer

    4.3. Name: root
         Description: the root identification number (default to 1 for the primary root of a plant)
         Type: integer

    4.4. Name: order
         Description: the root branching order
         Type: integer

    4.5. Name: parentroot
         Description: the identification number of the parent root
         Type: integer

    4.6. Name: time
         Description: the root system age
         Unit: days
         Type: integer

    4.7. Name: DBase
         Description: the distance between the branching point to the parent root base
         Unit: cm
         Type: numeric

    4.8. Name: length
         Description: the root length
         Unit: cm
         Type: numeric

    4.9. Name: mean.diameter
         Description: the average root diameter
         Unit: cm
         Type: numeric

    4.10. Name: sd.diameter
          Description: the standard deviation of the root diameter
          Unit: cm
          Type: numeric

    4.11. Name: nlat
          Description: the number of lateral roots
          Type: integer

    4.12. Name: branching.angle
          Description: the root branching angle
          Unit: degrees
          Type: numeric

    4.13. Name: tortuosity
          Description: the tortuosity (i.e., the ratio of the root length to the Euclidean distance between the branching point and the apex of the root)
          Unit: -
          Type: numeric

    4.14. Name: surface
          Description: the root surface area
          Unit: cm^2
          Type: numeric

    4.15. Name: volume
          Description: the root volume
          Unit: cm^3
          Type: numeric

    4.16. Name: lauz
          Description: the length of the unbranched apical zone
          Unit: cm
          Type: numeric

5. Missing data codes:
        NA        Not available