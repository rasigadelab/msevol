#include "event_loop.h"

/******************************************************************/
/* MAIN EVENT LOOP */

#include <iostream>

/* Run all events once */
void EventLoop::step() {

	/* TODO randomize event order too ?*/

	cellBirthAll();
	cellDeathAll();
	plasmidTransferAll();
	plasmidLossAll();
	cellDiffusionAll();

	u.shuffle();
}

/* Run all events */
void EventLoop::run(const size_t duration) {
	/* Register barcodes of all elements before startup */
	//u.register_all_barcodes();
	u.mset_register_all();

	for (size_t i = 0; i < duration; i++) {

		//cout << ".";

		//cout << endl << "Iter = " << i << endl;
		//for (auto it : u.vertices<Gene>()) {
		//	cout << "Gene " << it.second->id << " count = " << it.second->count() << endl;
		//}
		//for (auto it : u.vertices<Plasmid>()) {
		//	cout << "Plasmid " << it.second->id << " count = " << it.second->count() << endl; 
		//}
		//for (auto it : u.vertices<Chromosome>()) {
		//	cout << "Chromosome " << it.second->id << " count = " << it.second->count() << endl;
		//}
		//for (auto& it : u.vertices<Cell>()) {
		//	cout << "Cell " << it.id << " count = " << it.data.count() << endl;
		//}
		//for (auto it : u.vertices<Patch>()) {
		//	cout << "Patch " << it.second->id << " count = " << it.second->count() << endl;
		//}
		//cout << endl;

		//cout << endl << "Iter = " << i << endl;
		//auto s = u.size();
		//cout << "|V| = " << s.first << " --- |E| = " << s.second << endl;

		step();
	}
}