# Animated Cartograms
R script to create animated cartograms
<br>
<br>

----- ENGLISH VERSION ----- <br>
<i>This work was carried out by Hugo Thomas during his training period supervised by Florent Demoraes, at the CNRS ESO research unit in Rennes (France), in September and October 2021.</i>

## Uncovering urban circadian pulses based on an animated cartogram: the example of Bogotá

This blog presents the evolution of the spatial distribution of an urban population over 24 hours. More precisely, the objective is to display on a dynamic map the population rebalancing that takes place in the city during the day, taking as a starting point the population distribution in the early morning before the citizens move. The example developed focuses on Bogotá, the capital of Colombia, gathering 9 million inhabitants. The cartographic technique used is that of an animated cartogram, which is still not widely used in geography. The data used are micro-data from the 2018 census (population at place of residence) and from the latest OD survey conducted on a sample of 20,000 households in 2019 (population at place of activity). This work thus presents the originality of combining socio-demographic and daily mobility data in a dynamic mode. The programming was implemented in the R environment and relies on several packages including <a href="https://cran.r-project.org/web/packages/cartogramR/" target="_new" rel="noopener"><i>cartogramR</i></a> developed by <a href="https://perso.univ-rennes2.fr/pierre-andre.cornillon" target="_new" rel="noopener">Pierre-André Cornillon</a> (UMR CNRS 6625 - IRMAR : Institut de Recherche Mathématique de Rennes - Université Rennes 2) and <a href="https://perso.univ-rennes2.fr/florent.demoraes" target="_new" rel="noopener">Florent Demoraes</a> (UMR CNRS 6590 - ESO : Espaces et Sociétés - Université Rennes 2). In a way of open science and reproducibility, the R markdown script is provided below. 

This study is part of the activities of the <a href="https://modural.hypotheses.org/le-projet" target="_new" rel="noopener">Modural program</a> funded by the French National Research Agency and is conducted in collaboration with the French National Institute of Demographic Studies (INED) and the French National Scientific Research Center (CNRS), which are currently working on the adaptation of the <a href="https://mobiliscope.cnrs.fr/en" target="_new" rel="noopener">Mobiliscope project</a> to Latin American metropolises. The results clearly show the persistence of the segregated socio-residential pattern of the city and the maintenance of a very marked spatial mismatch between places of residence and places of work, with the center polarizing most of the flows from the densely populated deprived outskirts.

#### Key words 
<i>Circadian urban pulses, daily population rebalancing, diurnal densities, spatial mismatch, micro data, animated cartogram, R script, Bogotá</i>

--> Access to the <a href="https://github.com/ESO-Rennes/Animated-Cartograms/blob/main/pulsations.Rmd" target="_new" rel="noopener"><strong>R markdown script</strong></a>

Acknowledgements: <i>Pierre-André Cornillon, Guillaume Le Roux, Aurélie Douet</i>

<br>
<br>
----- VERSION FRANCAISE ----- <br>

<i>Ce travail a été réalisé par Hugo Thomas lors de son stage encadré par Florent Demoraes, à l'unité de recherche CNRS ESO de Rennes (France), en septembre et octobre 2021.</i>

## Révéler les pulsations urbaines circadiennes à l’aide d’un cartogramme animé : l’exemple de Bogotá

Ce billet présente l’évolution de la distribution spatiale d’une population urbaine sur 24 heures. Plus précisément, l’objectif est de restituer sur une carte dynamique les rééquilibrages de population qui s’opèrent dans la ville pendant la journée en ayant comme point de départ la répartition de la population au petit matin avant que les citadins ne bougent. L’exemple développé porte sur Bogotá, capitale de la Colombie, rassemblant 9 millions d’habitants. La technique cartographique employée est celle d’un cartogramme animé, dont l’usage est encore peu répandu en géographie. Les données utilisées sont les micro-données du recensement de 2018 (population sur le lieu de résidence) et celles de la dernière enquête OD menée en 2019 auprès d’un échantillon de 20 000 ménages (population sur le lieu d’activité). Ce travail présente ainsi l’originalité de combiner en mode dynamique des données socio-démographiques et de mobilité quotidienne. La programmation a été mise en œuvre dans l’environnement R et fait appel à plusieurs packages dont <a href="https://cran.r-project.org/web/packages/cartogramR/" target="_new" rel="noopener"><i>cartogramR</i></a> développé par <a href="https://perso.univ-rennes2.fr/pierre-andre.cornillon" target="_new" rel="noopener">Pierre-André Cornillon</a> (UMR CNRS 6625 - IRMAR : Institut de Recherche Mathématique de Rennes - Université Rennes 2) et <a href="https://perso.univ-rennes2.fr/florent.demoraes" target="_new" rel="noopener">Florent Demoraes</a>. Dans une perspective de science ouverte et de reproductibilité, le script R (markdown) est fourni plus bas.

Cette étude s’inscrit dans les activités du <a href="https://modural.hypotheses.org/le-projet" target="_new" rel="noopener">programme Modural</a> financé par l’Agence Nationale de la Recherche française et est menée en lien avec l’Institut National d’Études Démographiques de France (Ined) et le Centre National de la Recherche Scientifique de France (CNRS) qui travaillent actuellement sur l’extension du <a href="https://mobiliscope.cnrs.fr/fr" target="_new" rel="noopener">projet Mobiliscope</a> aux métropoles d’Amérique latine. Les résultats mettent bien en évidence la persistance du modèle socio-résidentiel ségrégatif de la ville et le maintien d’un découplage spatial très marqué entre lieux de résidence et lieux d’activité, le centre polarisant l’essentiel des flux en provenance de la périphérie populaire densément peuplée.

#### Mots-clefs
<i>Pulsations urbaines circadiennes, rééquilibrage populationnel journalier, densités diurnes, découplage spatial, micro-données, cartogramme animé, script R, Bogotá</i>

--> Accéder au <a href="https://github.com/ESO-Rennes/Animated-Cartograms/blob/main/pulsations.Rmd" target="_new" rel="noopener"><strong> fichier R markdown</strong></a>

Remerciements : <i>Pierre-André Cornillon, Guillaume Le Roux, Aurélie Douet</i>
