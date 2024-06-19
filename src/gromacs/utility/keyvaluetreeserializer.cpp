/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "gromacs/utility/keyvaluetreeserializer.h"

#include <cstdint>

#include <map>
#include <mutex>
#include <string>
#include <typeindex>
#include <utility>
#include <vector>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/iserializer.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

namespace gmx
{

namespace
{

class ValueSerializer
{
public:
    static void initSerializers();

    static void              serialize(const KeyValueTreeValue& value, ISerializer* serializer);
    static KeyValueTreeValue deserialize(ISerializer* serializer);

private:
    ValueSerializer() = delete;

    typedef void (*SerializerFunction)(const KeyValueTreeValue& value, ISerializer* serializer);
    typedef void (*DeserializerFunction)(KeyValueTreeValueBuilder* builder, ISerializer* serializer);

    struct Serializer
    {
        unsigned char        typeTag;
        SerializerFunction   serialize;
        DeserializerFunction deserialize;
    };
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    static std::mutex s_initMutex;
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    static std::map<std::type_index, Serializer> s_serializers;
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
    static std::map<unsigned char, DeserializerFunction> s_deserializers;
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::mutex ValueSerializer::s_initMutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::map<std::type_index, ValueSerializer::Serializer> ValueSerializer::s_serializers;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
std::map<unsigned char, ValueSerializer::DeserializerFunction> ValueSerializer::s_deserializers;

template<typename T>
struct SerializationTraits
{
};

template<>
struct SerializationTraits<KeyValueTreeObject>
{
    static void serialize(const KeyValueTreeObject& value, ISerializer* serializer)
    {
        int count = value.properties().size();
        serializer->doInt(&count);
        for (const auto& prop : value.properties())
        {
            serializer->doString(const_cast<std::string*>(&prop.key()));
            ValueSerializer::serialize(prop.value(), serializer);
        }
    }
    static void deserialize(KeyValueTreeValueBuilder* value, ISerializer* serializer)
    {
        KeyValueTreeObjectBuilder builder(value->createObject());
        deserializeObject(&builder, serializer);
    }
    static void deserializeObject(KeyValueTreeObjectBuilder* builder, ISerializer* serializer)
    {
        int         count = 0;
        std::string key;
        serializer->doInt(&count);
        for (int i = 0; i < count; ++i)
        {
            serializer->doString(&key);
            builder->addRawValue(key, ValueSerializer::deserialize(serializer));
        }
    }
};

template<>
struct SerializationTraits<KeyValueTreeArray>
{
    static void serialize(const KeyValueTreeArray& array, ISerializer* serializer)
    {
        int count = array.values().size();
        serializer->doInt(&count);
        for (const auto& value : array.values())
        {
            ValueSerializer::serialize(value, serializer);
        }
    }
    static void deserialize(KeyValueTreeValueBuilder* value, ISerializer* serializer)
    {
        KeyValueTreeArrayBuilder builder(value->createArray());
        int                      count = 0;
        serializer->doInt(&count);
        for (int i = 0; i < count; ++i)
        {
            builder.addRawValue(ValueSerializer::deserialize(serializer));
        }
    }
};

template<>
struct SerializationTraits<std::string>
{
    static void serialize(const std::string& value, ISerializer* serializer)
    {
        serializer->doString(const_cast<std::string*>(&value));
    }
    static void deserialize(KeyValueTreeValueBuilder* builder, ISerializer* serializer)
    {
        std::string value;
        serializer->doString(&value);
        builder->setValue<std::string>(value);
    }
};

template<>
struct SerializationTraits<bool>
{
    static void serialize(bool value, ISerializer* serializer) { serializer->doBool(&value); }
    static void deserialize(KeyValueTreeValueBuilder* builder, ISerializer* serializer)
    {
        bool value = false;
        serializer->doBool(&value);
        builder->setValue<bool>(value);
    }
};

template<>
struct SerializationTraits<char>
{
    static void serialize(char value, ISerializer* serializer) { serializer->doChar(&value); }
    static void deserialize(KeyValueTreeValueBuilder* builder, ISerializer* serializer)
    {
        char value = ' ';
        serializer->doChar(&value);
        builder->setValue<char>(value);
    }
};

template<>
struct SerializationTraits<unsigned char>
{
    static void serialize(unsigned char value, ISerializer* serializer)
    {
        serializer->doUChar(&value);
    }
    static void deserialize(KeyValueTreeValueBuilder* builder, ISerializer* serializer)
    {
        unsigned char value = '\0';
        serializer->doUChar(&value);
        builder->setValue<unsigned char>(value);
    }
};

template<>
struct SerializationTraits<int>
{
    static void serialize(int value, ISerializer* serializer) { serializer->doInt(&value); }
    static void deserialize(KeyValueTreeValueBuilder* builder, ISerializer* serializer)
    {
        int value = 0;
        serializer->doInt(&value);
        builder->setValue<int>(value);
    }
};

template<>
struct SerializationTraits<int64_t>
{
    static void serialize(int64_t value, ISerializer* serializer) { serializer->doInt64(&value); }
    static void deserialize(KeyValueTreeValueBuilder* builder, ISerializer* serializer)
    {
        int64_t value = 0;
        serializer->doInt64(&value);
        builder->setValue<int64_t>(value);
    }
};

template<>
struct SerializationTraits<float>
{
    static void serialize(float value, ISerializer* serializer) { serializer->doFloat(&value); }
    static void deserialize(KeyValueTreeValueBuilder* builder, ISerializer* serializer)
    {
        float value = 0;
        serializer->doFloat(&value);
        builder->setValue<float>(value);
    }
};

template<>
struct SerializationTraits<double>
{
    static void serialize(double value, ISerializer* serializer) { serializer->doDouble(&value); }
    static void deserialize(KeyValueTreeValueBuilder* builder, ISerializer* serializer)
    {
        double value = 0;
        serializer->doDouble(&value);
        builder->setValue<double>(value);
    }
};

//! Helper function for serializing values of a certain type.
template<typename T>
void serializeValueType(const KeyValueTreeValue& value, ISerializer* serializer)
{
    SerializationTraits<T>::serialize(value.cast<T>(), serializer);
}

#define SERIALIZER(tag, type)                                                       \
    {                                                                               \
        std::type_index(typeid(type)),                                              \
        {                                                                           \
            tag, &serializeValueType<type>, &SerializationTraits<type>::deserialize \
        }                                                                           \
    }

// static
void ValueSerializer::initSerializers()
{
    std::lock_guard<std::mutex> lock(s_initMutex);
    if (!s_serializers.empty())
    {
        return;
    }
    s_serializers = {
        SERIALIZER('O', KeyValueTreeObject),
        SERIALIZER('A', KeyValueTreeArray),
        SERIALIZER('s', std::string),
        SERIALIZER('b', bool),
        SERIALIZER('c', char),
        SERIALIZER('u', unsigned char),
        SERIALIZER('i', int),
        SERIALIZER('l', int64_t),
        SERIALIZER('f', float),
        SERIALIZER('d', double),
    };
    for (const auto& item : s_serializers)
    {
        s_deserializers[item.second.typeTag] = item.second.deserialize;
    }
}

void ValueSerializer::serialize(const KeyValueTreeValue& value, ISerializer* serializer)
{
    auto iter = s_serializers.find(value.type());
    GMX_RELEASE_ASSERT(iter != s_serializers.end(), "Unknown value type for serializization");
    unsigned char typeTag = iter->second.typeTag;
    serializer->doUChar(&typeTag);
    iter->second.serialize(value, serializer);
}

KeyValueTreeValue ValueSerializer::deserialize(ISerializer* serializer)
{
    unsigned char typeTag = 0;
    serializer->doUChar(&typeTag);
    auto iter = s_deserializers.find(typeTag);
    GMX_RELEASE_ASSERT(iter != s_deserializers.end(), "Unknown type tag for deserializization");
    KeyValueTreeValueBuilder builder;
    iter->second(&builder, serializer);
    return builder.build();
}

} // namespace

//! \cond libapi
void serializeKeyValueTree(const KeyValueTreeObject& root, ISerializer* serializer)
{
    GMX_RELEASE_ASSERT(!serializer->reading(), "Incorrect serializer direction");
    ValueSerializer::initSerializers();
    SerializationTraits<KeyValueTreeObject>::serialize(root, serializer);
}

KeyValueTreeObject deserializeKeyValueTree(ISerializer* serializer)
{
    GMX_RELEASE_ASSERT(serializer->reading(), "Incorrect serializer direction");
    ValueSerializer::initSerializers();
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder obj(builder.rootObject());
    SerializationTraits<KeyValueTreeObject>::deserializeObject(&obj, serializer);
    return builder.build();
}
//! \endcond

} // namespace gmx
