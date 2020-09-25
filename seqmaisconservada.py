


from dnaTrie.principal import trieNo


class indels:

    root = trieNo('+')

    especiesList = root.leEspecies()
    qtdeEspecies = root.contaEspecies()
    porceCon = root.porcConservacao(qtdeEspecies)
    locaisConservados = root.locaisConserv(porceCon)
    seqAminos = root.leSequencias()
    indices = root.agrupaSitesCons(locaisConservados)
    conservados = root.geraConserv(indices, seqAminos)
    predominantes = root.geraPred(conservados)
    #maisCons = root.maiorSequenciaConservada(conservados, seqAminos)




    # gera arquivo de saida

    arquivo = "saidaMaisCons - Tipo"+str(root.tipo)+".txt"
    with open(arquivo, "w") as arq:


        text09 = '\nSequência Mais Conservada\n\n\n'
        arq.write(text09)



        listaSaida = []
        for i  in range(len(predominantes)):
            for j in range(len(predominantes[i])):
                listaSaida.append(predominantes[i][j])
                arq.write(predominantes[i][j])


        arq.write('\n\n')
        arq.write('--------------------------------------------------------------------------------------------------------------------------')
        arq.write('\n\n')
        arq.close()
























