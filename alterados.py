


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




    # gera arquivo de saida

    arquivo = "saidaAlter - Tipo"+str(root.tipo)+".txt"
    with open(arquivo, "w") as arq:


        text09 = '\nSequências Conservadas e Alteradas\n'

        arq.write(text09)


        for j in range(len(indices)) :
            arq.write('\nSites:\n')
            arq.write(str(indices[j]))
            arq.write('\n')
            arq.write('\nSequência conservada mais comum:\n')
            arq.write(predominantes[j])
            arq.write('\n')





            arq.write('\nResumo do alinhamento na sequência:\n')
            listConser = root.compara(predominantes[j], conservados[j], especiesList)[4]
            arq.write('\n - Espécies Conservadas:\n')
            for item in range(len(listConser)):
                cons = str(listConser[item])
                arq.write(cons)
                arq.write('\n')

            valSitAlt = root.compara(predominantes[j], conservados[j], especiesList)[3]

            arq.write('\n - Espécies com Alterações:\n')
            arq.write('\nTotal:\n')
            arq.write(str(valSitAlt))
            arq.write('\n\n')

            list1 = root.compara(predominantes[j], conservados[j], especiesList)[0]
            list2 = root.compara(predominantes[j], conservados[j], especiesList)[1]
            list3 = root.compara(predominantes[j], conservados[j], especiesList)[2]

            root.ordena(list1, list2, list3)


            for item in range(len(list1)):

                seq =  str(list1[item])
                esp = str(list2[item])
                alt = str(list3[item])


                arq.write('Sequência: {} - Espécie: {} - Alterações: {}'.format(seq, esp, alt))

                arq.write('\n')

            arq.write('\n')


            arq.write('--------------------------------------------------------------------------------------------------------------------------')
            arq.write('\n\n')



    arq.close()
























