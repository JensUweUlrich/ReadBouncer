//Not run
// defeniere eine globale Variable int active_threads = 0 ; 
// Die Funktion Threading All 
threading_all{
    vec shapes = shape_generator();

    for(shape s in shapes )
    {
        while (active_threads == 25)
        
            sleep(1000); // damit ist 1 Sekunbde gemeint , nachschauen  --> In Header hinzufügen , nachschauen 

        
        thread t = thread(run_program ,s)
        active_threads++;
        t.start();
    }
}
// die Idee ist es , nach jedem Thread eine Puase zu haben von 1 Sekunde #+


// run_program soll wie folgt verändert werden : 
run_prpgram(shape s)
{
    CustomBloomFilter bf = CreatBloomFilter(s,...);
    .....
}
CreatBloomFilter(s,...)
{
    ComputeMinimizer(s,.....);
}
// TODO :
/*
1) add sleep to header ThreadPool // done
2) Creat BloomFilter soll keine Ausgaben mehr haben // nachher
3) ComputeMinimizer auch 
4) Defeniere Vector von Shapes als glbale Variable , dann die Shapes zu dem Hinzufügen dann verwenden // done
in Sketching und Minimizer 
5) Addiere eine neue Parameter zu run_programm, ComputeMinimizer, Sketching
*/ 
