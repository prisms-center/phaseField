/*
 * parallelNucleationList.h
 *
 *  Created on: Mar 10, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_PARALLELNUCLEATIONLIST_H_
#define INCLUDE_PARALLELNUCLEATIONLIST_H_

# include "nucleus.h"

template <int dim>
class parallelNucleationList
{
public:
	parallelNucleationList(std::vector<nucleus<dim> > _newnuclei);
	std::vector<nucleus<dim> > buildGlobalNucleiList(double min_dist_between_nuclei, double min_dist_between_OP, unsigned int old_num_nuclei);
	std::vector<nucleus<dim> > removeSubsetOfNuclei(std::vector<unsigned int> nuclei_to_remove, unsigned int nuclei_size);

protected:
	void sendUpdate (int procno) const;
	void receiveUpdate (int procno);
	void broadcastUpdate (int broadcastProc, int thisProc);
	void resolveNucleationConflicts (double min_dist_between_nuclei, double min_dist_between_OP, unsigned int old_num_nuclei);
	std::vector<nucleus<dim> > newnuclei;

};

#endif /* INCLUDE_PARALLELNUCLEATIONLIST_H_ */
