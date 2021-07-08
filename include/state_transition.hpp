//
// Created by L.Nagy on 15/05/2021.
//

#ifndef SD_COOLING_INCLUDE_STATE_TRANSITION_HPP
#define SD_COOLING_INCLUDE_STATE_TRANSITION_HPP

#include <cmath>
#include <stdexcept>
#include <utility>
#include <functional>

#include "formatter.hpp"

/**
 * This function is called when there are fewer new LEM states than old LEM
 * states.
 * @tparam StateList a list of states
 * @tparam OccupancyList
 * @tparam DistanceFunction
 * @param lems
 * @param occupancy
 * @param new_lems
 * @return
 */
template<typename State, typename Real>
std::pair<std::vector<State>, std::vector<Real>>
more_old_lems_than_new(
        const std::vector<State> &old_lems,
        const std::vector<Real> &old_occupancy,
        const std::vector<State> &new_lems
) {
    // Verify that old_lems and old_occupancy have the same size.
    if (old_lems.size() != old_occupancy.size()) {
        throw std::runtime_error(
                Formatter() << "Expected: old_lems and old_occupancy to be of same size."
        );
    }

    // Verify that there are more old states than new - this means we need
    // somehow to deal with the old states.
    if (new_lems.size() >= old_lems.size()) {
        throw std::runtime_error(
                Formatter() << "Expected fewer new LEM states than old"
        );
    }

    std::vector<State> lems;
    std::vector<Real> occupancy;

    // Find the new LEMs that correspond with the old LEMs and build the occupancy
    // vector in the correct manner.
    for (int new_idx = 0; new_idx < new_lems.size(); ++new_idx) {
        auto new_lem = new_lems[new_idx];
        // For the given new_lem find the corresponding index of the old_lem.
        int corresponding_old_lem_idx = -1;
        Real distance = 1E100;
        for (int old_idx = 0; old_idx < old_lems.size(); ++old_idx) {
            auto old_lem = old_lems[old_idx];
            Real new_distance = std::abs(old_lem - new_lem);
            if (new_distance < distance) {
                distance = new_distance;
                corresponding_old_lem_idx = old_idx;
            }
        }
        if (corresponding_old_lem_idx >= 0) {
            lems.push_back(new_lem);
            occupancy.push_back(old_occupancy[corresponding_old_lem_idx]);
        }
    }

    // Since there are fewer LEMs now, we must make sure that the occupancy vector
    // sums to unity.
    // TODO: this is normalized based on a very simple scheme, perhaps something
    //       more sophisticated is required here.
    Real sum = 0;
    for (auto rho : occupancy) {
        sum += rho;
    }
    for (int i = 0; i < occupancy.size(); ++i) {
        occupancy[i] /= sum;
    }

    return {lems, occupancy};
}

template<typename State, typename Real>
std::pair<std::vector<State>, std::vector<Real>>
more_new_lems_than_old(
        const std::vector<State> &old_lems,
        const std::vector<Real> &old_occupancy,
        const std::vector<State> &new_lems
) {
    // Verify that old_lems and old_occupancy have the same size.
    if (old_lems.size() != old_occupancy.size()) {
        throw std::runtime_error(
                Formatter() << "Expected: old_lems and old_occupancy to be of same size."
        );
    }

    // Verify that there are more or equal new states than old.
    if (new_lems.size() < old_lems.size()) {
        throw std::runtime_error(
                Formatter() << "Expected fewer new LEM states than old"
        );
    }

    std::vector<State> lems;
    std::vector<Real> occupancy;

    // Find new LEMs that correspond with the old LEMs and build the occupancy
    // in the correct way.
    std::vector<bool> new_lem_idx_picked(new_lems.size(), false);
    for (int old_idx = 0; old_idx < old_lems.size(); ++old_idx) {
        auto old_lem = old_lems[old_idx];
        // For the given old_lem find the corresponding index of the new_lem.
        int corresponding_new_lem_idx = -1;
        Real distance = 1E100;
        for (int new_idx = 0; new_idx < new_lems.size(); ++new_idx) {
            // Only process new_lems that remain.
            if (!new_lem_idx_picked[new_idx]) {
                auto new_lem = new_lems[new_idx];
                Real new_distance = std::abs(old_lem - new_lem);
                if (new_distance < distance) {
                    distance = new_distance;
                    corresponding_new_lem_idx = new_idx;
                }
            }
        }
        if (corresponding_new_lem_idx >= 0) {
            new_lem_idx_picked[corresponding_new_lem_idx] = true;
            lems.push_back(new_lems[corresponding_new_lem_idx]);
            occupancy.push_back(old_occupancy[old_idx]);
        }
    }

    // Now if there remain any unprocessed new_lems, set their occupancy to zero.
    for (int new_idx = 0; new_idx < new_lems.size(); ++new_idx) {
        if (!new_lem_idx_picked[new_idx]) {
            lems.push_back(new_lems[new_idx]);
            occupancy.push_back(0);
        }
    }

    return {lems, occupancy};
}

template<typename State, typename Real>
std::pair<std::vector<State>, std::vector<Real>>
same_no_of_new_lems_as_old(
        const std::vector<State> &old_lems,
        const std::vector<Real> &old_occupancy,
        const std::vector<State> &new_lems
) {
    // Verify that old_lems and old_occupancy have the same size.
    if (old_lems.size() != old_occupancy.size()) {
        throw std::runtime_error(
                Formatter() << "Expected: old_lems and old_occupancy to be of same size."
        );
    }

    // Verify that the number of old LEM states and new LEM states are the same.
    if (new_lems.size() != old_lems.size()) {
        throw std::runtime_error(
                Formatter() << "Expected fewer new LEM states than old"
        );
    }

    std::vector<State> lems;
    std::vector<Real> occupancy;

    // Find new LEMs that correspond with the old LEMs and build the occupancy
    // in the correct way.
    std::vector<bool> new_lem_idx_picked(new_lems.size(), false);
    for (int old_idx = 0; old_idx < old_lems.size(); ++old_idx) {
        auto old_lem = old_lems[old_idx];
        // For the given old_lem find the corresponding index of the new_lem.
        int corresponding_new_lem_idx = -1;
        Real distance = 1E100;
        for (int new_idx = 0; new_idx < new_lems.size(); ++new_idx) {
            // Only process new_lems that remain.
            if (!new_lem_idx_picked[new_idx]) {
                auto new_lem = new_lems[new_idx];
                Real new_distance = std::abs(old_lem - new_lem);
                if (new_distance < distance) {
                    distance = new_distance;
                    corresponding_new_lem_idx = new_idx;
                }
            }
        }
        if (corresponding_new_lem_idx >= 0) {
            new_lem_idx_picked[corresponding_new_lem_idx] = true;
            lems.push_back(new_lems[corresponding_new_lem_idx]);
            occupancy.push_back(old_occupancy[old_idx]);
        }
    }

    return {lems, occupancy};
}

#endif //SD_COOLING_INCLUDE_STATE_TRANSITION_HPP
