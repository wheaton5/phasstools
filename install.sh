
cd purge_dups/src
make
cd ../..
cd FASTK
make
cd ../het_snp_kmers
cargo build --release
cd ../molecule_kmers
cargo build --release
cd ../phasst_phase
cargo build --release
cd ../phasst_scaff
cargo build --release
cd ..

