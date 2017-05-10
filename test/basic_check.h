#include "sstream"
#include "string"
#include "cmath"

#include "basic_types.h"
#include "inputinstance.h"

typedef enum
{
    SUCCESS,
    FAILURE,
    NO_ZERO
} ReturnType;

struct ReturnMessage
{
    ReturnMessage(const ReturnType &t, const std::string &m)
        : type(t)
        , message(m)
    {
    }

    ReturnMessage(const ReturnType &t)
        : type(t)
        , message("")
    {
    }

    ReturnMessage()
        : type(ReturnType::SUCCESS)
        , message("")
    {
    }

    ReturnMessage(const ReturnMessage& other)
        : type(other.type)
        , message(other.message)
    {
    }

    ReturnType type;
    std::string message;
};


const InputInstance makeAllEqualInstance(const unsigned int num_chr, const unsigned int num_sam,
                                         const unsigned int num_seg, const double value);
const InputInstance makeCompleteIntInstance(const unsigned int num_chr, const unsigned int num_sam);
const InputInstance makeCompleteFracInstance(const unsigned int num_chr, const unsigned int num_sam);

