/***************************************************************************
*                                                                          *
*   Copyright (C) 2018 by David B. Blumenthal                              *
*                                                                          *
*   This file is part of GEDLIB.                                           *
*                                                                          *
*   GEDLIB is free software: you can redistribute it and/or modify it      *
*   under the terms of the GNU Lesser General Public License as published  *
*   by the Free Software Foundation, either version 3 of the License, or   *
*   (at your option) any later version.                                    *
*                                                                          *
*   GEDLIB is distributed in the hope that it will be useful,              *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           *
*   GNU Lesser General Public License for more details.                    *
*                                                                          *
*   You should have received a copy of the GNU Lesser General Public       *
*   License along with GEDLIB. If not, see <http://www.gnu.org/licenses/>. *
*                                                                          *
***************************************************************************/

/*!
 * @file compression.ipp
 * @brief ged::COMPRESSION class definition.
 */

#ifndef SRC_EDIT_COSTS_COMPRESSION_IPP_
#define SRC_EDIT_COSTS_COMPRESSION_IPP_

namespace ged {

template<>
COMPRESSION<GXLLabel, GXLLabel>::
~COMPRESSION() {}

template<>
COMPRESSION<GXLLabel, GXLLabel>::
COMPRESSION(double node_ins_cost, double node_del_cost, double node_rel_cost, double edge_ins_cost, double edge_del_cost, double edge_rel_cost,  double edge_rel_cost_id):
node_ins_cost_{node_ins_cost},
node_del_cost_{node_del_cost},
node_rel_cost_{node_rel_cost},
edge_ins_cost_{edge_ins_cost},
edge_del_cost_{edge_del_cost},
edge_rel_cost_{edge_rel_cost},
edge_rel_cost_id_{edge_rel_cost_id} {}

template<>
double
COMPRESSION<GXLLabel, GXLLabel>::
node_ins_cost_fun(const GXLLabel & node_label) const {
	return node_ins_cost_;
}

template<>
double
COMPRESSION<GXLLabel, GXLLabel>::
node_del_cost_fun(const GXLLabel & node_label) const {
	return node_del_cost_;
}

template<>
double
COMPRESSION<GXLLabel, GXLLabel>::
node_rel_cost_fun(const GXLLabel & node_label_1, const GXLLabel & node_label_2) const {
	bool is_zero=true;
	if (node_label_1.size()!=node_label_2.size()){
		is_zero=false;
	}
	else{
		for(auto const& entry : node_label_1){
			if(node_label_2.count(entry.first)<1){
				is_zero=false;
			}
			if(node_label_2.at(entry.first) != entry.second){
				is_zero=false;
			}
		}	
	}
	
	if (is_zero) {
		return 0.0;
	}
	return node_rel_cost_;
}

template<>
GXLLabel
COMPRESSION<GXLLabel, GXLLabel>::
median_node_label(const std::vector<GXLLabel> & node_labels) const {
	// Construct histogram.
	std::map<GXLLabel, std::size_t> hist;
	for (const auto & label : node_labels) {
		if (hist.find(label) == hist.end()) {
			hist[label] = 1;
		}
		else {
			hist[label]++;
		}
	}

	// Return the label that appears most frequently.
	std::size_t best_count{0};
	GXLLabel median_label;
	for (const auto & label_count : hist) {
		if (label_count.second > best_count) {
			best_count = label_count.second;
			median_label = label_count.first;
		}
	}
	return median_label;
}

template<>
double
COMPRESSION<GXLLabel, GXLLabel>::
edge_ins_cost_fun(const GXLLabel & edge_label) const {
	return edge_ins_cost_;
}

template<>
double
COMPRESSION<GXLLabel, GXLLabel>::
edge_del_cost_fun(const GXLLabel & edge_label) const {
	// Need to differenciate edge removal because of node removal or without node removal
	// This however will be done before calling this function, as the edge label is not enough to differentiate the two cases 
	return edge_del_cost_;
}

template<>
double
COMPRESSION<GXLLabel, GXLLabel>::
edge_rel_cost_fun(const GXLLabel & edge_label_1, const GXLLabel & edge_label_2) const {
	bool is_id=true;
	if (edge_label_1.size()!=edge_label_2.size()){
		is_id=false;
	}
	else{
		for(auto const& entry : edge_label_1){
			if(edge_label_2.count(entry.first)<1){
				is_id=false;
			}
			if(edge_label_2.at(entry.first) != entry.second){
				is_id=false;
			}
		}
	}
	
	if (is_id) {
		return edge_rel_cost_id_;
	}
	return edge_rel_cost_;
}

template<>
GXLLabel
COMPRESSION<GXLLabel, GXLLabel>::
median_edge_label(const std::vector<GXLLabel> & edge_labels) const {
	// Construct histogram.
	std::map<GXLLabel, std::size_t> hist;
	for (const auto & label : edge_labels) {
		if (hist.find(label) == hist.end()) {
			hist[label] = 1;
		}
		else {
			hist[label]++;
		}
	}

	// Return the label that appears most frequently.
	std::size_t best_count{0};
	GXLLabel median_label;
	for (const auto & label_count : hist) {
		if (label_count.second > best_count) {
			best_count = label_count.second;
			median_label = label_count.first;
		}
	}
	return median_label;
}

}

#endif /* SRC_EDIT_COSTS_COMPRESSION_IPP_ */
