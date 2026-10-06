#!/usr/bin/env python3

"""MIT License

Copyright (c) 2025 James Eaton, Andrew Baldwin (University of Oxford)
              2025, Bind Research

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE."""

"""
Working out which maxima down the bore dimension belong to which peak of a 2D
reference plane.

Peaks which are close to one another in the plane share the same region of the
3D data, so a maximum found down the bore of one of them can belong to its
neighbour instead. The functions here take the maxima found down the bore of
each reference peak and work out which of them each peak most likely owns,
saying when the answer is not clear enough to rely on.

Everything here works on numbers alone so that it can be checked without the
graphical interface.
"""

import numpy as np


# How a peak came by the position it has down the bore dimension
CLEAR = "clear"
AMBIGUOUS = "ambiguous"
RESOLVED = "resolved"

# How sure the choice between two peaks has to be before it is treated as
# clear. A candidate which both peaks have an equal claim to scores 0.5, so a
# margin of 0.1 treats anything from 0.4 to 0.6 as ambiguous
DEFAULT_MARGIN = 0.1

# How close two peaks have to be in the plane, in points, before a maximum
# down the bore of one of them could belong to the other
DEFAULT_CLOSENESS_POINTS = 4


def find_default_closeness(ppms_0, ppms_1, points=DEFAULT_CLOSENESS_POINTS):
    """
    How close two peaks have to be in the plane before their maxima down the
    bore can be confused, as a distance in the units of each plane axis. The
    distance is a few points of the data, which follows the spacing of the
    spectrum rather than assuming a number of ppm.
    """
    distances = []

    for ppms in [ppms_0, ppms_1]:
        try:
            spacing = abs(float(ppms[1]) - float(ppms[0]))
        except (IndexError, TypeError, ValueError):
            spacing = 0.0

        distances.append(spacing * points)

    return distances


def find_close_peaks(positions, distances):
    """
    Which peaks are close enough to one another in the plane that a maximum
    down the bore of one of them could belong to the other.

    positions is a list of (shift1, shift2) and distances is how close the
    peaks have to be in each of those two dimensions. The answer is a list
    holding, for each peak, the positions in the list of its neighbours.
    """
    neighbours = [[] for position in positions]

    for i, first in enumerate(positions):
        for j, second in enumerate(positions):
            if i == j:
                continue

            close = True
            for dimension in [0, 1]:
                distance = abs(float(first[dimension]) - float(second[dimension]))
                if distance > abs(float(distances[dimension])):
                    close = False
                    break

            if close == True:
                neighbours[i].append(j)

    return neighbours


def find_column_fit(shifts, low, high) -> float:
    """
    How much of a column of chemical shifts falls inside the range of an axis,
    as a fraction of the peaks. A column which belongs to that axis fits nearly
    all of the way.
    """
    shifts = [shift for shift in shifts]

    if len(shifts) == 0:
        return 0.0

    low, high = min(float(low), float(high)), max(float(low), float(high))

    inside = 0
    for shift in shifts:
        try:
            if low <= float(shift) <= high:
                inside += 1
        except (TypeError, ValueError):
            continue

    return inside / len(shifts)


def find_axis_assignment(columns, ranges, minimum=0.5, score=None, margin=0.05):
    """
    Which column of a peaklist belongs to which axis of the spectrum being
    shown.

    columns holds the three columns of chemical shifts as they are in the
    peaklist, and ranges the low and high chemical shift of each axis, in the
    order the peaklist should end up in. The answer says which column to take
    for each axis, so a peaklist picked from a different projection is plotted
    the right way round.

    Two dimensions of a spectrum can cover the same range of chemical shifts,
    as the two nitrogens of an (H)N(CA)NH do, and then the ranges alone cannot
    say which column belongs to which of them. Where more than one arrangement
    fits the ranges as well as the others, score is used to choose between
    them: it is given an arrangement and says how well the peaks land on the
    data with the columns that way round, which the ranges cannot know.

    None is returned when no arrangement of the columns fits the spectrum, as
    happens when the peaklist belongs to another spectrum altogether.
    """
    if len(columns) != 3 or len(ranges) != 3:
        return None

    orders = [
        (0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0),
    ]

    fitting = []

    for order in orders:
        fits = [
            find_column_fit(columns[order[axis]], ranges[axis][0], ranges[axis][1])
            for axis in range(3)
        ]

        # More than half of each column has to fall within its axis, otherwise
        # the peaklist does not belong to this spectrum. A peaklist with a few
        # peaks outside the range is still used, as a peaklist can hold peaks
        # which the projection being shown does not cover
        if min(fits) <= minimum:
            continue

        fitting.append((list(order), sum(fits)))

    if len(fitting) == 0:
        return None

    best_fit = max([total for order, total in fitting])

    # The arrangements which fit the ranges about as well as one another
    close = [order for order, total in fitting if best_fit - total <= margin]

    if len(close) == 1 or score == None:
        # The first of them keeps the peaklist as it is where there is nothing
        # to choose between them
        return close[0]

    best = None
    best_score = None

    for order in close:
        try:
            value = float(score(order))
        except (TypeError, ValueError):
            continue

        if best_score == None or value > best_score:
            best_score = value
            best = order

    if best == None:
        return close[0]

    return best


def find_trace_maxima(trace, minimum=0.0):
    """
    Every maximum of a trace down the bore dimension: each point which is
    larger than the points on either side of it and at least as large as the
    minimum given.

    The peak picker only reports the maxima which clear the threshold it is
    given, so a bore which holds more than one resonance can have maxima which
    it does not report. Those maxima are worth knowing about, as one of them can
    be the one a peak really belongs to.
    """
    indexes = []

    values = [abs(float(value)) for value in trace]

    for index in range(len(values)):
        if values[index] < abs(float(minimum)):
            continue

        before = values[index - 1] if index > 0 else None
        after = values[index + 1] if index < len(values) - 1 else None

        if before != None and values[index] < before:
            continue
        if after != None and values[index] <= after:
            continue

        indexes.append(index)

    return indexes


def find_extra_candidates(trace, shifts, minimum, known, separation=1):
    """
    The maxima of a trace down the bore which the peak picker did not report,
    which are the other positions a peak could have. A maximum next to one which
    is already known is left out, as it is the same resonance.
    """
    extra = []

    for index in find_trace_maxima(trace, minimum):
        if any(abs(index - held) <= separation for held in known):
            continue
        if any(abs(index - held["bore_index"]) <= separation for held in extra):
            continue

        try:
            shift = float(shifts[index])
        except (IndexError, TypeError, ValueError):
            continue

        extra.append(
            {
                "bore_index": int(index),
                "shift": shift,
                "intensity": float(trace[index]),
                "extra": True,
            }
        )

    return extra


def find_projection_fit(projection, intensity):
    """
    How well the intensity of a maximum down the bore fits the strength of a
    peak in the projection.

    The projection holds the largest intensity of each position down the bore,
    so a maximum belonging to a peak is about as strong as that peak is in the
    projection. A maximum which is much stronger than a peak is in the
    projection more likely belongs to a stronger neighbour, and one which is
    much weaker more likely belongs to a weaker one.
    """
    projection = abs(float(projection))
    intensity = abs(float(intensity))

    if projection == 0 or intensity == 0:
        return 0.0

    return min(projection, intensity) / max(projection, intensity)


def find_intensity(data, plane, bore_index):
    """
    The intensity of the 3D data at a position in the plane and a position
    down the bore dimension. The data is held as [bore][plane][plane].
    """
    try:
        return float(data[int(bore_index)][int(plane[0])][int(plane[1])])
    except (IndexError, TypeError, ValueError):
        return 0.0


def find_confidence(
    data,
    plane,
    neighbour_planes,
    bore_index,
    projection=None,
    neighbour_projections=None,
    intensity=None,
):
    """
    How strongly a maximum down the bore belongs to a peak rather than to one
    of its neighbours in the plane.

    Two things are weighed up. The intensity of the 3D data at the peak is
    compared with the intensity of the same slice of the bore at each of its
    neighbours, and, when the strengths of the peaks in the projection are
    known, how well the maximum fits the strength of this peak in the
    projection is compared with how well it fits its neighbours. A peak which
    holds all of the claim scores 1, one which shares it equally with a
    neighbour scores 0.5, and one whose neighbour has the better claim scores
    below 0.5.
    """
    own = abs(find_intensity(data, plane, bore_index))

    others = [
        abs(find_intensity(data, neighbour, bore_index))
        for neighbour in neighbour_planes
    ]

    if len(others) == 0:
        # Nothing else is close enough to have a claim on this maximum
        return 1.0

    best_other = max(others)

    if own + best_other == 0:
        intensity_score = 0.0
    else:
        intensity_score = own / (own + best_other)

    if projection == None or neighbour_projections == None or intensity == None:
        return intensity_score

    if len(neighbour_projections) == 0:
        return intensity_score

    # How well the maximum fits the strength of each peak in the projection
    own_fit = find_projection_fit(projection, intensity)
    other_fit = max(
        [
            find_projection_fit(neighbour, intensity)
            for neighbour in neighbour_projections
        ]
    )

    if own_fit + other_fit == 0:
        return intensity_score

    projection_score = own_fit / (own_fit + other_fit)

    # The two are given equal weight: where the data at the peak and the
    # strength of the peak in the projection agree the answer is clear, and
    # where they disagree the maximum is worth looking at
    return (intensity_score + projection_score) / 2


def find_candidate_scores(data, peaks, neighbours, margin=DEFAULT_MARGIN):
    """
    Score every maximum found down the bore of every peak.

    peaks is a list of dictionaries holding:
      plane      - the position of the peak in the data as (index, index)
      candidates - the maxima found down its bore, each holding "bore_index",
                   "shift" and "intensity"

    Each candidate gains a "confidence" saying how strongly it belongs to this
    peak rather than to a neighbour, and "ambiguous" saying whether the two
    have too similar a claim on it to be sure.
    """
    scored = []

    for i, peak in enumerate(peaks):
        neighbour_planes = [peaks[j]["plane"] for j in neighbours[i]]
        neighbour_projections = [peaks[j].get("projection") for j in neighbours[i]]
        if any(value == None for value in neighbour_projections):
            neighbour_projections = None

        candidates = []
        for candidate in peak["candidates"]:
            confidence = find_confidence(
                data,
                peak["plane"],
                neighbour_planes,
                candidate["bore_index"],
                peak.get("projection"),
                neighbour_projections,
                candidate.get("intensity"),
            )

            held = dict(candidate)
            held["confidence"] = confidence
            held["ambiguous"] = abs(confidence - 0.5) <= margin
            candidates.append(held)

        scored.append(candidates)

    return scored


def find_assignments(
    data, peaks, neighbours, expected=0, margin=DEFAULT_MARGIN
):
    """
    Work out which maxima down the bore belong to each peak of the reference
    plane.

    The maxima of each peak are scored against its neighbours, the ones the
    peak most likely owns are kept, and the rest are held as the other
    positions the peak could have. Where expected is more than zero that many
    maxima are kept for each peak, and a peak with fewer than that is noted:
    a peak can be missing a maximum for good reasons, such as being next to a
    proline or at the end of the chain.

    The answer holds, for each peak, a list of the maxima it has been given
    (each with its position, intensity, confidence and whether it is
    ambiguous), the other positions it could have had, and a note of anything
    the user should look at.
    """
    scored = find_candidate_scores(data, peaks, neighbours, margin)

    assignments = []

    for i, candidates in enumerate(scored):
        # The maxima this peak has the strongest claim to come first
        ordered = sorted(
            candidates,
            key=lambda candidate: (
                candidate["confidence"],
                abs(candidate["intensity"]),
            ),
            reverse=True,
        )

        # A maximum which a neighbour holds more of belongs to the neighbour
        claimed = [
            candidate for candidate in ordered if candidate["confidence"] >= 0.5
        ]
        given_away = [
            candidate for candidate in ordered if candidate["confidence"] < 0.5
        ]

        if expected > 0:
            kept = claimed[:expected]
        else:
            kept = claimed

        others = claimed[len(kept):] + given_away

        notes = []
        if any(candidate["ambiguous"] == True for candidate in kept):
            notes.append(AMBIGUOUS)
        if expected > 0 and len(kept) < expected:
            notes.append("{} of {} found".format(len(kept), expected))

        assignments.append(
            {
                "kept": kept,
                "others": others,
                "notes": notes,
                "candidates": ordered,
            }
        )

    return assignments


def find_ambiguity_text(notes) -> str:
    """
    How the notes of a peak are shown in the peaklist table and written into a
    peaklist file.
    """
    if len(notes) == 0:
        return ""

    # Only the first letter is changed, so that the names of peaks inside a note
    # keep the case they are written in
    return ", ".join(
        [str(note)[:1].upper() + str(note)[1:] for note in notes]
    )


def find_alternatives_text(others) -> str:
    """
    How the other positions a peak could have had are shown in the peaklist
    table and written into a peaklist file. They are separated by semicolons so
    that they stay in one column.
    """
    if len(others) == 0:
        return ""

    return ";".join(["{:.5f}".format(candidate["shift"]) for candidate in others])


def read_alternatives(text) -> list:
    """
    The other positions a peak could have had, read back from a peaklist file.
    """
    alternatives = []

    for part in str(text).split(";"):
        part = part.strip()
        if part == "" or part == "nan":
            continue
        try:
            alternatives.append(float(part))
        except ValueError:
            continue

    return alternatives


def find_plane_distance(first, second, distances):
    """
    How far apart two positions are in the plane, counted in how much of the
    distance allowed in each dimension is used up. A position which uses up all
    of one dimension is a distance of 1 away, so anything above 1 is too far.
    """
    used = []

    for dimension in [0, 1]:
        allowed = abs(float(distances[dimension]))
        gap = abs(float(first[dimension]) - float(second[dimension]))

        if allowed == 0:
            used.append(0.0 if gap == 0 else float("inf"))
        else:
            used.append(gap / allowed)

    return max(used), (used[0] ** 2 + used[1] ** 2) ** 0.5


def find_template_matches(picked, template, distances):
    """
    Relate peaks picked in the whole of a 3D to the peaks of a template
    peaklist, using the two dimensions they share: the plane.

    picked holds the peaks which were found, each with "shift1", "shift2",
    "shift3" and "intensity"; template holds the peaks to relate them to, each
    with "name", "shift1" and "shift2"; and distances says how far apart the
    shared dimensions can be for the peaks to be the same one.

    The answer says, for each picked peak, which template peaks are close
    enough to it, nearest first. A picked peak which is close to more than one
    of them cannot be told apart by the plane alone, and one which is close to
    none of them is not in the template at all.
    """
    matches = []

    for peak in picked:
        position = (peak["shift1"], peak["shift2"])

        owners = []
        for index, held in enumerate(template):
            inside, distance = find_plane_distance(
                position, (held["shift1"], held["shift2"]), distances
            )
            if inside <= 1:
                owners.append((distance, index))

        owners.sort()

        matches.append(
            {
                "owners": [index for distance, index in owners],
                "distances": [distance for distance, index in owners],
            }
        )

    return matches


def find_template_assignments(picked, template, distances, expected=0):
    """
    Work out which of the peaks picked in a 3D belongs to each peak of a
    template peaklist.

    Each template peak is given the peaks which sit on it in the plane, the
    nearest first, keeping as many as are expected where a number is given. A
    picked peak which more than one template peak could claim is noted, along
    with the names of the others it could belong to, and the peaks which no
    template peak is near are gathered separately so that they can be looked at
    too.
    """
    matches = find_template_matches(picked, template, distances)

    assignments = [
        {"kept": [], "others": [], "notes": [], "shared": []}
        for held in template
    ]
    unmatched = []

    for i, match in enumerate(matches):
        if len(match["owners"]) == 0:
            unmatched.append(i)
            continue

        for place, owner in enumerate(match["owners"]):
            shared = [
                template[other]["name"]
                for other in match["owners"]
                if other != owner
            ]

            assignments[owner]["kept"].append(
                {
                    "picked": i,
                    "distance": match["distances"][place],
                    "shift1": picked[i]["shift1"],
                    "shift2": picked[i]["shift2"],
                    "shift": picked[i]["shift3"],
                    "intensity": picked[i]["intensity"],
                    "ambiguous": len(shared) > 0,
                    "shared": shared,
                }
            )

    for i, assignment in enumerate(assignments):
        # The peaks which sit closest to the template peak come first, and the
        # strongest of those at the same distance
        assignment["kept"].sort(
            key=lambda peak: (peak["distance"], -abs(peak["intensity"]))
        )

        if expected > 0 and len(assignment["kept"]) > expected:
            assignment["others"] = assignment["kept"][expected:]
            assignment["kept"] = assignment["kept"][:expected]

        notes = []
        if any(peak["ambiguous"] == True for peak in assignment["kept"]):
            notes.append(AMBIGUOUS)
        if expected > 0 and len(assignment["kept"]) < expected:
            notes.append("{} of {} found".format(len(assignment["kept"]), expected))
        if len(assignment["kept"]) == 0 and expected == 0:
            notes.append("none found")

        assignment["notes"] = notes

    return assignments, unmatched
