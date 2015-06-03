"""
Provides classes that take protocol requests, send that request to
the server, and write a particular genomics file type with the results.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import collections
import string

import pysam

import ga4gh.datamodel.reads as reads


class AbstractConverter(object):
    """
    Abstract base class for converter classes
    """
    def __init__(self, httpClient):
        self._httpClient = httpClient


##############################################################################
# SAM
##############################################################################


class SamException(Exception):
    """
    Something that went wrong during converting a SAM file
    """


class SamConverter(AbstractConverter):
    """
    Converts a request to a SAM file
    """
    def __init__(self, httpClient, searchReadsRequest, outputFile,
                 binaryOutput):
        super(SamConverter, self).__init__(httpClient)
        self._searchReadsRequest = searchReadsRequest
        self._outputFile = outputFile
        self._binaryOutput = binaryOutput

    def convert(self):
        header = self._getHeader()
        targetIds = self._getTargetIds(header)
        # pysam can't write to file streams (except for stdout)
        # http://pysam.readthedocs.org/en/latest/usage.html#using-streams
        if self._binaryOutput:
            flags = "wb"
        else:
            flags = "wh"  # h for header
        fileString = "-"
        if self._outputFile is not None:
            fileString = self._outputFile
        alignmentFile = pysam.AlignmentFile(
            fileString, flags, header=header)
        iterator = self._httpClient.searchReads(self._searchReadsRequest)
        for read in iterator:
            alignedSegment = SamLine.toAlignedSegment(read, targetIds)
            alignmentFile.write(alignedSegment)
        alignmentFile.close()

    def _getHeader(self):
        # TODO where to get actual values for header?
        # need some kind of getReadGroup(readGroupId) method in protocol
        # just add these dummy lines for now
        header = {
            'HD': {'VN': '1.0'},
            'SQ': [
                {'LN': 1575, 'SN': 'chr1'},
                {'LN': 1584, 'SN': 'chr2'},
            ],
        }
        return header

    def _getTargetIds(self, header):
        # this seems to be how pysam sets the target ids
        targetIds = collections.defaultdict(int)
        targetId = 0
        if 'SQ' in header:
            headerLines = header['SQ']
            for headerLine in headerLines:
                refName = headerLine['SN']
                targetIds[refName] = targetId
                targetId += 1
        return targetIds


class SamLine(object):
    """
    Methods for processing a line in a SAM file
    """
    _encoding = 'utf8'

    # see tables in SAM spec, section 1.5
    _tagReservedFieldPrefixes = set(["X", "Y", "Z", ])
    _tagIntegerFields = set([
        "AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", "HI", "IH", "MQ",
        "NH", "NM", "OP", "PQ", "SM", "TC", "UQ", ])
    _tagStringFields = set([
        "BC", "BQ", "CC", "CO", "CQ", "CS", "CT", "E2", "FS", "LB", "MC",
        "MD", "OQ", "OC", "PG", "PT", "PU", "QT", "Q2", "R2", "RG", "RT",
        "SA", "U2", ])
    _tagIntegerArrayFields = set(["FZ", ])

    def __init__(self):
        raise SamException("SamLine can't be instantiated")

    @classmethod
    def toAlignedSegment(cls, read, targetIds):
        ret = pysam.AlignedSegment()
        # QNAME
        ret.query_name = read.fragmentName.encode(cls._encoding)
        # SEQ
        ret.query_sequence = read.alignedSequence.encode(cls._encoding)
        # FLAG
        ret.flag = cls.toSamFlag(read)
        # RNAME
        refName = read.alignment.position.base.referenceName
        ret.reference_id = targetIds[refName]
        # POS
        ret.reference_start = int(read.alignment.position.base.position)
        # MAPQ
        ret.mapping_quality = read.alignment.mappingQuality
        # CIGAR
        ret.cigar = cls.toCigar(read)
        # RNEXT
        nextRefName = read.nextMatePosition.base.referenceName
        ret.next_reference_id = targetIds[nextRefName]
        # PNEXT
        ret.next_reference_start = int(read.nextMatePosition.base.position)
        # TLEN
        ret.template_length = read.fragmentLength
        # QUAL
        ret.query_qualities = read.alignedQuality
        ret.tags = cls.toTags(read)
        return ret

    @classmethod
    def toSamFlag(cls, read):
        flag = 0
        if read.numberReads:
            reads.SamFlags.setFlag(
                flag, reads.SamFlags.NUMBER_READS)
        if read.properPlacement:
            reads.SamFlags.setFlag(
                flag, reads.SamFlags.PROPER_PLACEMENT)
        if read.readNumber:
            reads.SamFlags.setFlag(
                flag, reads.SamFlags.READ_NUMBER_ONE)
            reads.SamFlags.setFlag(
                flag, reads.SamFlags.READ_NUMBER_TWO)
        if read.secondaryAlignment:
            reads.SamFlags.setFlag(
                flag, reads.SamFlags.SECONDARY_ALIGNMENT)
        if read.failedVendorQualityChecks:
            reads.SamFlags.setFlag(
                flag, reads.SamFlags.FAILED_VENDOR_QUALITY_CHECKS)
        if read.duplicateFragment:
            reads.SamFlags.setFlag(
                flag, reads.SamFlags.DUPLICATE_FRAGMENT)
        if read.supplementaryAlignment:
            reads.SamFlags.setFlag(
                flag, reads.SamFlags.SUPPLEMENTARY_ALIGNMENT)
        return flag

    @classmethod
    def toCigar(cls, read):
        cigarTuples = []
        for gaCigarUnit in read.alignment.cigar:
            operation = reads.SamCigar.ga2int(gaCigarUnit.operation)
            length = int(gaCigarUnit.operationLength)
            cigarTuple = (operation, length)
            cigarTuples.append(cigarTuple)
        return tuple(cigarTuples)

    @classmethod
    def _parseTagValue(cls, tag, value):
        if tag[0] in cls._tagReservedFieldPrefixes:
            # user reserved fields... not really sure what to do here
            return value[0].encode(cls._encoding)
        elif tag in cls._tagIntegerFields:
            return int(value[0])
        elif tag in cls._tagStringFields:
            return value[0].encode(cls._encoding)
        elif tag in cls._tagIntegerArrayFields:
            return [int(integerString) for integerString in value]
        else:
            raise SamException("unrecognized tag '{}'".format(tag))

    @classmethod
    def toTags(cls, read):
        tags = []
        for tag, value in read.info.items():
            val = cls._parseTagValue(tag, value)
            tagTuple = (tag, val)
            tags.append(tagTuple)
        retval = tuple(tags)
        return retval


##############################################################################
# VCF
##############################################################################


class VcfException(Exception):
    pass


class VcfConverter(AbstractConverter):
    """
    Converts a request to a VCF file
    """
    def __init__(self, httpClient, outputStream, searchVariantSetsRequest,
                 searchVariantsRequest):
        raise NotImplementedError


##############################################################################
# Lastgraph
##############################################################################


class LastgraphConverter(AbstractConverter):
    """
    Converts a request to a lastgraph file
    """
    def __init__(self, httpClient, outputStream, searchSequencesRequest,
        searchJoinsRequest, saveSequence=True):
        """
        Given the HTTP client used to make the requests, a stream to write the
        lastgraph-format output to, a request for all the sequences in a
        ReferenceSet, and a request for all the joins in a ReferenceSet, creates
        a converter that converts the graph defined by the sequences and joins
        to lastgraph format. Every sequence will become a node, and every join
        an arc.
        
        If saveSequence is false, don't download the sequence data.
        """

        super(LastgraphConverter, self).__init__(httpClient)

        # Save the arguments
        self._outputStream = outputStream
        self._searchSequencesRequest = searchSequencesRequest
        self._searchJoinsRequest = searchJoinsRequest
        self._saveSequence = saveSequence
        
        # We need a cache of sequence bases for populating our nodes, so we
        # don't make multiple requests for split nodes. Holds full sequences by
        # ID.
        self._sequence_cache = {}
        
        # We need this for taking reverse complements. Handle all the ambiguity
        # codes as defined at <http://www.boekhoff.info/?pid=data&dat=fasta-
        # codes>
        self._complement_table = string.maketrans(
            "ACGTNKMRYSWBVHDacgtnkmryswbvhd", "TGCANMKYRSWVBDHtgcanmkyrswvbdh")
        
    def _get_sequence(self, sequence_id):
        """
        Return the bases for the sequence with the given ID. Downloads it if it
        has not been downloaded already.
        """
        
        if not self._sequence_cache.has_key(sequence_id):
            # We need to do the download, and convert to bytes when we get it.
            self._sequence_cache[sequence_id] = \
                self._httpClient.getSequenceBases(sequence_id).sequence.encode(
                "utf8")
               
        # Give back what we have  
        return self._sequence_cache[sequence_id]
        
    def _reverse_complement(self, string):
        """
        Return the reverse complement of a DNA string.
        
        """
        
        # Complement and reverse it.
        return str.translate(string, self._complement_table)[::-1]
        
        
    def convert(self):
        """
        Actually do the conversion.
        """

        # Get an iterator over the sequences
        sequences = self._httpClient.searchSequences(
            self._searchSequencesRequest)

        # Get an iterator over the joins
        joins = self._httpClient.searchJoins(self._searchJoinsRequest)
        
        # TODO: Download all the sequences and joins and convert from string
        # graph to side graph (split nodes on joins).
        
        # Keep a set of (base, face) pairs on which joins are incedent in every
        # sequence, in a dict by sequence ID.
        sequence_breakers = collections.defaultdict(set)
        
        # Keep the downloaded joins so we don't have to get them again. Holds
        # pairs of (sequence, (index, is_right_side)) tuples. TODO: memory
        # usage?
        downloaded_joins = []
        
        def to_endpoint(side):
            """
            Turn a join Side into a sequence ID and an (index, is_right_side)
            pair.
            
            """
            
            return side.base.sequenceId, (int(side.base.position), 
                side.strand == "NEG_STRAND")
                
        def implicit_partner(endpoint):
            """
            Given a (index, is_right_side) pair, produce the implicitly joined
            partner.
            
            """
            
            if endpoint[1]:
                # Right sides need left sides
                return (endpoint[0] + 1, False)
            else:
                # Left sides need right sides
                return (endpoint[0] - 1, True)
                
        
        # Fill them in
        for join in joins:
            # For each join
            
            # Save it
            downloaded_joins.append((to_endpoint(join.side1),
                to_endpoint(join.side2)))
            
            # Break sequences at both ends, by adding (index, is_right_side)
            # pairs to the set of necessary breakpoints. We do this format so
            # they sort correctly.
            for seq, endpoint in [to_endpoint(side) 
                for side in [join.side1, join.side2]]:
                
                # For each place on each sequence the join touches, break the
                # sequence there.
                sequence_breakers[seq].add(endpoint)
                
                # Make sure to add the partner endpoint for the implicit join we
                # are making. It may be off the sequence; we will deal with that
                # later.
                sequence_breakers[seq].add(implicit_partner(endpoint))
        
        # We need to know how long each node is. This stores (sequence_id,
        # start, length) tuples for each node. The node number is 1 + the index
        # here.
        nodes = []
        
        # What's the total length of all nodes?
        total_length = 0
        
        # We now need to build this join (sequence_id, (index, is_right_side))
        # tuple to node number (node index + 1, negated for right side) dict.
        node_for_endpoint = {}
        
        for sequence in sequences:
            # For every sequence object, we're going to split it (now that we
            # have the length).
            
            # Add the left end and right end endpoints.
            sequence_breakers[sequence.id].add((0, False))
            sequence_breakers[sequence.id].add((int(sequence.length) - 1, True))
            
            # Grab all the places we need to have ends after we break the
            # sequence. We are guaranteed proper segments (at least one) with
            # endpoints on each end that are within the bounds of the sequence.
            endpoints = sorted((endpoint 
                for endpoint in sequence_breakers[sequence.id] 
                if endpoint[0] >= 0 and endpoint[0] < int(sequence.length)))
            
            for i in xrange(0, len(endpoints) - 1):
                if endpoints[i][1]:
                    # We're looking at an implicit join, since we have a right
                    # and then (presumably) a left. Add it.
                    downloaded_joins.append(((sequence.id, endpoints[i]),
                        (sequence.id, endpoints[i + 1])))
                else:
                    # We're looking at a segment
                    
                    # Where does this segment/node start in the sequence, and
                    # how long does it run?
                    start = endpoints[i][0]
                    length = endpoints[i + 1][0] - start + 1
                    
                    # Note the contribution to the total length
                    total_length += length
                    
                    # Add the node
                    nodes.append((sequence.id, start, length))
                    
                    # Point to it for when we do the joins. Left side gets node
                    # number (index + 1), right side gets node number negated.
                    # When we output the join, we will have to additionally
                    # negate its source.
                    node_for_endpoint[(sequence.id, endpoints[i])] = \
                        len(nodes)
                    node_for_endpoint[(sequence.id, endpoints[i + 1])] = \
                        -len(nodes)
                
        # Now we have the whole string graph in memory and can spit it out.
                
        # Put out a lastgraph header.
        self._outputStream.write("{}\t0\t0\t0\n".format(len(nodes)))

        for index, (sequence_id, start, length) in enumerate(nodes):
            # Announce every node with its number (= index + 1) and length
            
            self._outputStream.write("NODE\t{}\t{}\t0\n".format(index + 1,
                length))
            
            if self._saveSequence:
            
                # Get the node sequence by slicing the (cached) full sequence
                node_sequence = self._get_sequence(
                    sequence_id)[start:start + length]
                
                # Save the forward sequence
                self._outputStream.write(node_sequence)
                self._outputStream.write("\n")
                
                # Save the reverse sequence
                self._outputStream.write(self._reverse_complement(node_sequence))
                self._outputStream.write("\n")
                
            else:
                # Make up sequences
                
                for i in xrange(2):
                    for j in xrange(length):
                        self._outputStream.write("N")
                    self._outputStream.write("\n")
                
        # Now the arcs.
        for end1, end2 in downloaded_joins:
            # For each join we recorded, make an arc between those pre-negated
            # node numbers. We have to additionally negate the first one again,
            # because positive numbers all around correspond to a right to left
            # join.
            self._outputStream.write("ARC\t{}\t{}\n".format(
                -node_for_endpoint[end1], node_for_endpoint[end2]))

