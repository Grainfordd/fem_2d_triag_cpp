import gmsh

def crear_malla():
    gmsh.initialize()

    gm = gmsh.model.geo

    alto = 4
    ancho = 4

    p1 =  gm.addPoint(0, 0, 0) 
    p2 =  gm.addPoint(ancho, 0, 0) 
    p3 =  gm.addPoint(ancho, alto, 0) 
    p4 =  gm.addPoint(0, alto, 0) 

    l1 = gm.add_line(p1, p2)
    l2 = gm.add_line(p2, p3)
    l3 = gm.add_line(p3, p4)
    l4 = gm.add_line(p4, p1)

    n = 10
    m = 10

    cl1 = gm.add_curve_loop([l1, l2, l3, l4])

    for tag in [l1, l3]:
        gm.mesh.setTransfiniteCurve(tag, n)

    for tag in [l2, l4]:
        gm.mesh.setTransfiniteCurve(tag, m)


    s1 = gm.addPlaneSurface([cl1])
    gm.mesh.setTransfiniteSurface(s1)

    gmsh.model.addPhysicalGroup(2, [s1], 101) # Darle tag, para añadir physical name después
    gmsh.model.setPhysicalName(2, 101, 'Superficie')

    gm.synchronize()
    gmsh.model.mesh.generate(2)

    gmsh.write('malla_python.msh')
    info = gmsh.model.mesh.getNodes()

    # gmsh.fltk.run()

    return info, ancho

def crear_fuerzas(ancho, info):
    tags, coords, _ = info 
    valor_fuerza = 10e3

    nodos_extremo_derecho = []

    for i, tag in enumerate(tags):
        x = coords[i*3]
        y = coords[i*3 + 1]
        z = coords[i*3 + 2]

        tol = 1e-6
        if abs(x - ancho) < tol:
            nodos_extremo_derecho.append((tag, x, y, z))

    n = len(nodos_extremo_derecho)
    fuerzas_por_nodo = valor_fuerza/n


    with open('fuerzas.csv', 'w') as archivo:
        archivo.write('Nodo, dir, valor\n')

        for nodo_info in nodos_extremo_derecho:
            tag = nodo_info[0]
            archivo.write(f'{tag}, 0, {fuerzas_por_nodo}\n')

        print("Archivo 'fuerzas.csv' creado exitosamente")

def crear_fixed(info):
    node_tags, node_coords, _ = info
    
    # Encontrar nodos en el extremo izquierdo (x = 0)
    min_x = 0.0
    tolerancia = 1e-6
    
    nodos_extremo_izquierdo = []
    
    for i, tag in enumerate(node_tags):
        x = node_coords[i*3]
        y = node_coords[i*3 + 1]
        z = node_coords[i*3 + 2]
        
        # Verificar si el nodo está en el extremo izquierdo
        if abs(x - min_x) < tolerancia:
            nodos_extremo_izquierdo.append((tag, x, y, z))
    
    
    
    with open('disp.csv', 'w') as archivo:
        archivo.write('Nodo,dir,valor\n')
        
        for nodo_info in nodos_extremo_izquierdo:
            tag = nodo_info[0]
            # Aplicar condición de frontera en dirección x (dir = 0)
            archivo.write(f'{tag},0,0\n')
            # Aplicar condición de frontera en dirección y (dir = 1)
            archivo.write(f'{tag},1,0\n')
    
    print("Archivo 'disp.csv' creado exitosamente")
    




info, ancho = crear_malla()
print()
crear_fuerzas(ancho, info)
crear_fixed(info)
