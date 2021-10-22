# Animated Cartograms
R script to create animated cartograms<br>


https://user-images.githubusercontent.com/45881567/138525838-d1a5fc91-60fd-487c-940a-da44d83f53d9.mp4



----- ENGLISH VERSION ----- <br>
<i>This work was carried out by <a href="https://fr.linkedin.com/in/ht128" target="_blank">Hugo Thomas</a> during his training period supervised by <a href="https://perso.univ-rennes2.fr/florent.demoraes" target="_blank" rel="noopener">Florent Demoraes</a>, at the CNRS ESO research unit in Rennes (France), in September and October 2021.</i>

## Uncovering urban circadian pulses based on an animated cartogram: the example of Bogotá

This blog presents the evolution of the spatial distribution of an urban population over 24 hours. More precisely, the objective is to display on a dynamic map the population rebalancing that takes place in the city during the day, taking as a starting point the population distribution in the early morning before the citizens move. The example developed focuses on Bogotá, the capital of Colombia, gathering 9 million inhabitants. The cartographic technique used is that of an animated cartogram, which is still not widely used in geography. The data used are micro-data from the 2018 census (population at place of residence) and from the latest OD survey conducted on a sample of 20,000 households in 2019 (population at place of activity). This work thus presents the originality of combining socio-demographic and daily mobility data in a dynamic mode. The programming was implemented in the R environment and relies on several packages including <a href="https://cran.r-project.org/web/packages/cartogramR/" target="_blank" rel="noopener"><i>cartogramR</i></a> developed by <a href="https://perso.univ-rennes2.fr/pierre-andre.cornillon" target="_blank" rel="noopener">Pierre-André Cornillon</a> (UMR CNRS 6625 - IRMAR : Institut de Recherche Mathématique de Rennes - Université Rennes 2) and <a href="https://perso.univ-rennes2.fr/florent.demoraes" target="_blank" rel="noopener">Florent Demoraes</a> (UMR CNRS 6590 - ESO : Espaces et Sociétés - Université Rennes 2). In a way of open science and reproducibility, the R markdown script is provided below. 

This study is part of the activities of the <a href="https://modural.hypotheses.org/le-projet" target="_blank" rel="noopener">Modural program</a> funded by the French National Research Agency and is conducted in collaboration with the French National Institute of Demographic Studies (INED) and the French National Scientific Research Center (CNRS), which are currently working on the adaptation of the <a href="https://mobiliscope.cnrs.fr/en" target="_blank" rel="noopener">Mobiliscope project</a> to Latin American metropolises. The results clearly show the persistence of the segregated socio-residential pattern of the city and the maintenance of a very marked spatial mismatch between places of residence and places of work, with the center polarizing most of the flows from the densely populated deprived outskirts.

#### Key words 
<i>Circadian urban pulses, daily population rebalancing, diurnal densities, spatial mismatch, micro data, animated cartogram, R script, Bogotá</i>

--> Access to the <a href="https://github.com/ESO-Rennes/Animated-Cartograms/blob/main/pulsations.Rmd" target="_blank" rel="noopener"><strong>R markdown script</strong></a>

Acknowledgements: <i>Pierre-André Cornillon, Guillaume Le Roux, Aurélie Douet</i>

<br>
<br>
----- VERSION FRANÇAISE ----- <br>

<i>Ce travail a été réalisé par <a href="https://fr.linkedin.com/in/ht128" target="_blank" rel="noopener">Hugo Thomas</a> lors de son stage encadré par <a href="https://perso.univ-rennes2.fr/florent.demoraes" target="_blank" rel="noopener">Florent Demoraes</a>, à l'unité de recherche CNRS ESO de Rennes (France), en septembre et octobre 2021.</i>

## Révéler les pulsations urbaines circadiennes à l’aide d’un cartogramme animé : l’exemple de Bogotá

Ce billet présente l’évolution de la distribution spatiale d’une population urbaine sur 24 heures. Plus précisément, l’objectif est de restituer sur une carte dynamique les rééquilibrages de population qui s’opèrent dans la ville pendant la journée en ayant comme point de départ la répartition de la population au petit matin avant que les citadins ne bougent. L’exemple développé porte sur Bogotá, capitale de la Colombie, rassemblant 9 millions d’habitants. La technique cartographique employée est celle d’un cartogramme animé, dont l’usage est encore peu répandu en géographie. Les données utilisées sont les micro-données du recensement de 2018 (population sur le lieu de résidence) et celles de la dernière enquête OD menée en 2019 auprès d’un échantillon de 20 000 ménages (population sur le lieu d’activité). Ce travail présente ainsi l’originalité de combiner en mode dynamique des données socio-démographiques et de mobilité quotidienne. La programmation a été mise en œuvre dans l’environnement R et fait appel à plusieurs packages dont <a href="https://cran.r-project.org/web/packages/cartogramR/" target="_blank" rel="noopener"><i>cartogramR</i></a> développé par <a href="https://perso.univ-rennes2.fr/pierre-andre.cornillon" target="_blank" rel="noopener">Pierre-André Cornillon</a> (UMR CNRS 6625 - IRMAR : Institut de Recherche Mathématique de Rennes - Université Rennes 2) et <a href="https://perso.univ-rennes2.fr/florent.demoraes" target="_blank" rel="noopener">Florent Demoraes</a> (UMR CNRS 6590 - ESO : Espaces et Sociétés - Université Rennes 2). Dans une perspective de science ouverte et de reproductibilité, le script R (markdown) est fourni plus bas.

Cette étude s’inscrit dans les activités du <a href="https://modural.hypotheses.org/le-projet" target="_blank" rel="noopener">programme Modural</a> financé par l’Agence Nationale de la Recherche française et est menée en lien avec l’Institut National d’Études Démographiques de France (Ined) et le Centre National de la Recherche Scientifique de France (CNRS) qui travaillent actuellement sur l’extension du <a href="https://mobiliscope.cnrs.fr/fr" target="_blank" rel="noopener">projet Mobiliscope</a> aux métropoles d’Amérique latine. Les résultats mettent bien en évidence la persistance du modèle socio-résidentiel ségrégatif de la ville et le maintien d’un découplage spatial très marqué entre lieux de résidence et lieux d’activité, le centre polarisant l’essentiel des flux en provenance de la périphérie populaire densément peuplée.

#### Mots-clefs
<i>Pulsations urbaines circadiennes, rééquilibrage populationnel journalier, densités diurnes, découplage spatial, micro-données, cartogramme animé, script R, Bogotá</i>

--> Accéder au <a href="https://github.com/ESO-Rennes/Animated-Cartograms/blob/main/pulsations.Rmd" target="_blank" rel="noopener"><strong> fichier R markdown</strong></a>

Remerciements : <i>Pierre-André Cornillon, Guillaume Le Roux, Aurélie Douet</i><br>
<br>
<br>

----- VERSIÓN ESPAÑOLA ----- <br>
<i>Este trabajo fue realizado por <a href="https://fr.linkedin.com/in/ht128" target="_blank">Hugo Thomas</a> durante su práctica supervisada por <a href="http://www.ifea.org.pe/investigadores/florent-demoraes/" target="_blank" rel="noopener">Florent Demoraes</a>, en la unidad de investigación del CNRS ESO en Rennes (Francia), en septiembre y octubre de 2021.</i>
<br>
## Revelar los pulsos circadianos urbanos con un cartograma animado: el ejemplo de Bogotá

Este post presenta la evolución de la distribución espacial de una población urbana a lo largo de 24 horas. Más concretamente, el objetivo es reproducir en un mapa dinámico el reequilibrio de la población que tiene lugar en la ciudad durante el día, tomando como punto de partida la distribución de la población a primera hora de la mañana, antes de que los habitantes salgan de su casa. El ejemplo desarrollado se refiere a Bogotá, la capital de Colombia, con 9 millones de habitantes. La técnica cartográfica utilizada es la de un cartograma animado, que todavía no se utiliza mucho en geografía. Los datos utilizados son los microdatos del censo de 2018 (población en el lugar de residencia) y los de la última encuesta OD realizada en 2019 sobre una muestra de 20.000 hogares (población en el lugar de actividad). Este trabajo presenta así la originalidad de combinar datos sociodemográficos y de movilidad diaria en modo dinámico. La programación se implementó en el entorno R y utiliza varios paquetes, entre ellos <a href="https://cran.r-project.org/web/packages/cartogramR/" target="_blank" rel="noopener"><i>cartogramR</i></a> desarrollado por <a href="https://perso.univ-rennes2.fr/pierre-andre.cornillon" target="_blank" rel="noopener">Pierre-André Cornillon</a> (UMR CNRS 6625 - IRMAR : Institut de Recherche Mathématique de Rennes - Université Rennes 2) y <a href="http://www.ifea.org.pe/investigadores/florent-demoraes/" target="_blank" rel="noopener">Florent Demoraes</a> (UMR CNRS 6590 - ESO : Espaces et Sociétés - Université Rennes 2). En una perspectiva de ciencia abierta y reproducibilidad, el script R (markdown) se proporciona a continuación.

Este estudio forma parte de las actividades del <a href="https://modural.hypotheses.org/le-projet-modural/el-proyecto" target="_blank" rel="noopener">programa Modural</a>, financiado por la Agencia Nacional de Investigación de Francia, y se realiza en colaboración con el Instituto Nacional de Estudios Demográficos (INED) y el Centro Nacional de Investigación Científica (CNRS), que actualmente trabajan en la ampliación del <a href="https://mobiliscope.cnrs.fr/en" target="_blank" rel="noopener"> proyecto Mobiliscope </a> a las ciudades latinoamericanas. Los resultados muestran claramente la persistencia del modelo socio-residencial segregado de la ciudad y el mantenimiento de una disociación espacial muy marcada entre los lugares de residencia y los lugares de actividad, con el centro polarizando la mayor parte de los flujos procedentes de la periferia popular densamente poblada.

#### Palabras clave
<i>pulsos urbanos circadianos, reequilibrio diario de la población, densidades diurnas, disociación espacial, microdatos, cartograma animado, script R, Bogotá</i>

--> Acceder al <a href="https://github.com/ESO-Rennes/Animated-Cartograms/blob/main/pulsations.Rmd" target="_blank" rel="noopener"><strong>archivo R markdown</strong></a>

Agradecimientos : <i>Pierre-André Cornillon, Guillaume Le Roux, Aurélie Douet</i>





<br>
</br> <strong>References</strong>

Bertin J (1967). Sémiologie graphique: les diagrammes, les réseaux, les cartes. Mouton & Gauthier-Villars, Paris-La Haye.

Bertin J (1983). Semiology of graphics: diagrams, networks, maps. The University of Wisconsin Press, Madison, (trans. W. j. Berg). ISBN 0-299-09060-4.

Cauvin C, Escobar F, Serradj A (2007). Cartographie thématique 2 – Des transformations incontournables. Traité IGAT  – Information Géographique et Aménagement du Territoire ; Aspects fondamentaux de l’analyse spatiale. Hermès-Lavoisier.

Cauvin C, Escobar F, Serradj A (2010).  Cartography and the Impact of the Quantitative Revolution. John Wiley & Sons, Inc. https://onlinelibrary.wiley.com/doi/book/10.1002/9781118558126.

Demoraes F., Bouquet M., Mericskay B.,(2021) – How visually effective are animated cartograms? Potential improvements based on the example of segregation in Bogotá (1993-2005), M@ppemonde. DOI:10.4000/mappemonde.5928 - https://hal.archives-ouvertes.fr/hal-03152983 

Demoraes F., Bouquet M., Mericskay B.,(2021) – L’efficacité visuelle des cartogrammes animés en question - Une piste d’amélioration à travers l’exemple de la ségrégation à Bogotá (1993-2005), M@ppemonde. DOI:10.4000/mappemonde.5813 - https://hal.archives-ouvertes.fr/hal-03029241 

Dougenik JA, Chrisman NR, Niemeyer DR (1985). “An  Algorithm  To Construct Continuous Area Cartograms.”   The Professional Geographer, 37(1), 75–81. https://onlinelibrary.wiley.com/doi/abs/10.1111/j.0033-0124.1985.00075.x

Gastner MT,  Newman MEJ  (2004).   “Diffusion-Based Method  for  Producing Density-Equalizing Maps.”  Proceedings of the National Academy of Sciences, 101(20), 7499–7504. https://doi.org/10.1073/pnas.0400280101

Gastner MT,  Seguy V, More P (2018). “Fast Flow-Based Algorithm  for Creating Density-Equalizing Map Projections.”  Proceedings of the National Academy of Sciences, 115(10), E2156–E2164. https://doi.org/10.1073/pnas.1712674115

Nusrat S, Kobourov S (2016). “The State of the Art  in Cartograms.”  Computer Graphics Forum, 35(3), 619–642. https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.12932

Pebesma E (2018). “Simple Features for R: Standardized Support for Spatial Vector Data.” The R Journal, 10(1), 439–446. doi:10.32614/RJ-2018-009. https://journal.r-project.org/archive/2018/RJ-2018-009/index.html

Tobler W (2004). “Thirty Five Years of Computer Cartograms.”  Annals of the Association of American Geographers, 94(1), 58–73. https://onlinelibrary.wiley.com/doi/full/10.1111/j.1467-8306.2004.09401004.x
