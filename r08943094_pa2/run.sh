make clean; make
./bin/fp 0.5 ../input_pa2/xerox.block ../input_pa2/xerox.nets ../output_pa2/xerox.rpt -gui xerox.png | tee log/xerox.log
./bin/fp 0.5 ../input_pa2/hp.block ../input_pa2/hp.nets ../output_pa2/hp.rpt -gui hp.png | tee log/hp.log
./bin/fp 0.5 ../input_pa2/apte.block ../input_pa2/apte.nets ../output_pa2/apte.rpt -gui apte.png | tee log/apte.log
./bin/fp 0.5 ../input_pa2/ami33.block ../input_pa2/ami33.nets ../output_pa2/ami33.rpt -gui ami33.png | tee log/ami33.log
./bin/fp 0.5 ../input_pa2/ami49.block ../input_pa2/ami49.nets ../output_pa2/ami49.rpt -gui ami49.png | tee log/ami49.log
# bin/fp 0.5 ../input_pa2/test.block ../input_pa2/test.nets ../output_pa2/test.rpt -gui test.png