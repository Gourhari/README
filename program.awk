# 26 de julio de 2007
# Para correr el programa:
# $ awk -f script.awk <archivo_de_entrada> archivo_de_salida
# Archivo para extraer datos de energias, orbitales moleculares frontera
#proveniente de calculos Single Point
BEGIN{K=1}
{
  if ($1=="SCF" && $2=="Done:" ) 
  {
  energy=$5
  }

  while ( $1=="Alpha" && $2=="occ.") {
	  homo=$NF
          getline

	if ( $2=="virt."){
	j++
	lumo=$5
	    if (j==K){
	    printf "%-5s %-5s %-18s %+2s %-10s %+2s %-10s %s\n",
            j,"Energia=", energy, "E-HOMO=", homo, "E-LUMO=", lumo, K
	    K=K+1
            }
	}
  }
} 

