# Proteasome_APlab repo

Carpeta Heatmap de la Fig. 1:
 
-          Las muestras que hemos utilizado están en el archivo datasets.txt.
-          El punto de partida es Tex_EffectorMemory_CD8.rds, que contiene las T cells previamente clasificadas como Tex o EffectorMemory por ProjecTILs (merged Seurat objects). Os mandamos el script.R.
-          Inicialmente consideramos 45 genes del proteasoma, pero habría qué genes son los más informativos.
 
Carpeta Análisis UCell:
 
-          Intentamos refinar un poco la selección de muestras en los 4 tipos de cáncer inicial (Fig. 1).
-          También consideramos los datasets que nos sugeristeis de BRCA (discriminando entre HER2+, ER+ y TNBC) y melanoma (discriminando entre Anti-PD1, Anti-CTLA4 y Anti-PD1-CTLA4).
-          En el archivo proteasome.docx están las firmas de genes que finalmente consideramos en UCell.
 
Hipótesis: Nuestros análisis indican que, a nivel transcripcional, los genes que codifican para el proteasoma están sobreexpreados en Tex. Probablemente, la actividad proteasomal aumente como una estrategia para combatir el estrés proteotóxico que experimentan las células exhaustas. Nuestra hipótesis es que, aumentando la actividad del proteasoma en estados de diferenciación anteriores a Tex, podemos prevenir o retrasar la aparición del agotamiento en las células T.
 
Llegados a este punto, necesitamos vuestra experiencia para refinar los análisis y generar algunas figuras presentables en el paper.
 
Para ello nos gustaría que nos ayudarais a:
 
1)      Integrar los diferentes subtipos de CD8 (definidos por ProjecTILs) en un único análisis (pan cáncer). Eliminad o incluid las muestras, de estos u otros datasets, que consideréis oportuno.
2)      Definir qué genes son informativos para la firma del proteasoma y proyectar su expresión en el atlas del punto 1.
3)      Ver si tras el tratamiento con Anti-PD1 y/o Anti-CTLA4 (estudios de BRCA y melanoma) la expresión de la firma de genes (punto 2) disminuye.
 
