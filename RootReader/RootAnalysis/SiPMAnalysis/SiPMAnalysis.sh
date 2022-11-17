#!/bin/bash

g++ SiPMAnalysis.cpp `root-config --libs --cflags` -o sipmAnalysis

# Das ist nur eine Beispiel Datei, damit das geht, musst du die Datei von Julian mit dieser hier vergleichen und anpassen
# Gerade die Compilerzeile (g++ ...) ist vermutlich anders. Da muss hier noch was ergänzt werden um die Vasilov Funktion so benutzen zu können.
#
#Es ging mir darum, dir zu zeigen, wie man es ohne hardcoding macht, also der Teil hier unten sollte alle root Dateien im rootfiles Ordner einlesen und dann TEILWEISE parsen und ans Skript übergeben
#Das wird so oder ählich in fast allen von meinen Analyse Skripten gemacht, dort etwas zu suchen kann also nicht schaden.
#
# Ich konnte das alles nicht testen, weil nur wenig Zeit etc.
#
# Weißt du wie man debugged?
# Da das compiled wird, bleibt erstmal nichts anderes übrig als variablen in die Console zu loggen um zu überprüfen wie diese überhaupt aussehen.
# Hier in Bash wäre das sowieso wie :echo $line
#
# In C++ ist das auch sehr wichtig. Dort über cout<<"TEST"<<endl; etwas loggen. 
#
#
#Allgemein kannst du anhand des Skriptes von Julian gut üben/lernen wie man mit Problemen im Code umgeht. Schritt 1 wäre dabei immer das Verständnis. Du kannst nichts debuggen, was du nicht verstehst. 
#Der Großteil des Skripts sind irgendwelche Parametertests oder das übliche ROOT Plot Zeugs (was in jedem Skript so oder ähnlich auftaucht). Die Berechnungsroutine ist lächerlich wenig. Wenn du es verstanden hast, wärst du
#eigentlich locker in der Lage, diesen Teil zu kopieren und in ein neues Skript zu packen. 
# Alle Skripte machen mehr oder weniger das selbe Zeug, hast du einmal den Code verstanden, dann hast du quasi alles verstanden.



##################
### INITIALIZE ###
#################

  #Iteriere durch alle root Dateien im Ordner rootfiles -> ersetzt das nervige Hardcoden
 for file in "../rootfiles/"/*; do
	#without extension

	if [ ! -z "$1" ]; then

		if [[ $file == *"$1"* ]]; then
			echo $line
		else
			continue
		fi

	fi

	rootFilePath=$file

	prefix="../rootfiles//"
	runName=${file#$prefix} #Remove prefix
	suffix=".root"
	runName=${runName%$suffix} #Remove suffix

	runNr=$(echo "$runName" | sed -r 's/^([^.]+).*$/\1/; s/^[^0-9]*([0-9]+).*$/\1/')
	
  #Wenn was schief geht, schau dir an wie die Variable aussieht:
  echo $runNr $runName

	./sipmAnalysis $runName $runNr $rootFilePath 
done
 