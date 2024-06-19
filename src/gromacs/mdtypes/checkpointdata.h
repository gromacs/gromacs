/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
/*! \libinternal \file
 * \brief Provides the checkpoint data structure for the modular simulator
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdtypes
 */

#ifndef GMX_MODULARSIMULATOR_CHECKPOINTDATA_H
#define GMX_MODULARSIMULATOR_CHECKPOINTDATA_H

#include <cstdint>
#include <cstdio>

#include <optional>
#include <string>
#include <type_traits>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
class ISerializer;

/*! \libinternal
 * \brief The operations on CheckpointData
 *
 * This enum defines the two modes of operation on CheckpointData objects,
 * reading and writing. This allows to template all access functions, which
 * in turn enables clients to write a single function for read and write
 * access, eliminating the risk of having read and write functions getting
 * out of sync.
 *
 * \ingroup module_modularsimulator
 */
enum class CheckpointDataOperation
{
    Read,
    Write,
    Count
};

/*! \internal
 * \brief Get an ArrayRef whose const-ness is defined by the checkpointing operation
 *
 * \tparam operation  Whether we are reading or writing
 * \tparam T          The type of values stored in the ArrayRef
 * \param container   The container the ArrayRef is referencing to
 * \return            The ArrayRef
 *
 * \see ArrayRef
 *
 * \ingroup module_modularsimulator
 */
template<CheckpointDataOperation operation, typename T>
ArrayRef<std::conditional_t<operation == CheckpointDataOperation::Write || std::is_const<T>::value, const typename T::value_type, typename T::value_type>>
makeCheckpointArrayRef(T& container)
{
    return container;
}

/*! \internal
 * \brief Get an ArrayRef to a C array whose const-ness is defined by the checkpointing operation
 *
 * \tparam operation  Whether we are reading or writing
 * \tparam T          The type of values stored in the ArrayRef
 * \param begin       Pointer to the beginning of array.
 * \param size        Number of elements in array.
 * \return            The ArrayRef
 *
 * \see ArrayRef
 *
 * \ingroup module_modularsimulator
 */
template<CheckpointDataOperation operation, typename T>
ArrayRef<std::conditional_t<operation == CheckpointDataOperation::Write || std::is_const<T>::value, const T, T>>
makeCheckpointArrayRefFromArray(T* begin, size_t size)
{
    return ArrayRef<T>(begin, begin + size);
}

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Struct allowing to check if data is serializable through the KeyValueTree serializer
 *
 * This list of types is copied from ValueSerializer::initSerializers()
 * Having this here allows us to catch errors at compile time
 * instead of having cryptic runtime errors
 */
template<typename T>
struct IsSerializableType
{
    static bool const value = std::is_same<T, std::string>::value || std::is_same<T, bool>::value
                              || std::is_same<T, int>::value || std::is_same<T, int64_t>::value
                              || std::is_same<T, float>::value || std::is_same<T, double>::value;
};

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Struct allowing to check if enum has a serializable underlying type
 */
//! {
template<typename T, bool = std::is_enum<T>::value>
struct IsSerializableEnum
{
    static bool const value = IsSerializableType<std::underlying_type_t<T>>::value;
};
template<typename T>
struct IsSerializableEnum<T, false>
{
    static bool const value = false;
};
//! }

/*! \libinternal
 * \ingroup module_modularsimulator
 * \brief Data type hiding checkpoint implementation details
 *
 * This data type allows to separate the implementation details of the
 * checkpoint writing / reading from the implementation of the checkpoint
 * clients. Checkpoint clients interface via the methods of the CheckpointData
 * object, and do not need knowledge of data types used to store the data.
 *
 * Templating allows checkpoint clients to have symmetric (templated)
 * implementations for checkpoint reading and writing.
 *
 * CheckpointData objects are dispatched via [Write|Read]CheckpointDataHolder
 * objects, which interact with the checkpoint reading from / writing to
 * file.
 */

template<CheckpointDataOperation operation>
class CheckpointData;

//! Convenience shortcut for reading checkpoint data.
using ReadCheckpointData = CheckpointData<CheckpointDataOperation::Read>;
//! Convenience shortcut for writing checkpoint data.
using WriteCheckpointData = CheckpointData<CheckpointDataOperation::Write>;

template<>
class CheckpointData<CheckpointDataOperation::Read>
{
public:
    /*! \brief Read or write a single value from / to checkpoint
     *
     * Allowed scalar types include std::string, bool, int, int64_t,
     * float, double, or any enum with one of the previously mentioned
     * scalar types as underlying type. Type compatibility is checked
     * at compile time.
     *
     * \tparam operation  Whether we are reading or writing
     * \tparam T          The type of the value
     * \param key         The key to [read|write] the value [from|to]
     * \param value       The value to [read|write]
     */
    //! {
    template<typename T>
    std::enable_if_t<IsSerializableType<T>::value, void> scalar(const std::string& key, T* value) const;
    template<typename T>
    std::enable_if_t<IsSerializableEnum<T>::value, void> enumScalar(const std::string& key, T* value) const;
    //! }

    /*! \brief Read or write an ArrayRef from / to checkpoint
     *
     * Allowed types stored in the ArrayRef include std::string, bool, int,
     * int64_t, float, double, and gmx::RVec. Type compatibility is checked
     * at compile time.
     *
     * \tparam operation  Whether we are reading or writing
     * \tparam T          The type of values stored in the ArrayRef
     * \param key         The key to [read|write] the ArrayRef [from|to]
     * \param values      The ArrayRef to [read|write]
     */
    //! {
    // Read ArrayRef of scalar
    template<typename T>
    std::enable_if_t<IsSerializableType<T>::value, void> arrayRef(const std::string& key,
                                                                  ArrayRef<T>        values) const;
    // Read ArrayRef of RVec
    void arrayRef(const std::string& key, ArrayRef<RVec> values) const;
    //! }

    /*! \brief Read or write a tensor from / to checkpoint
     *
     * \tparam operation  Whether we are reading or writing
     * \param key         The key to [read|write] the tensor [from|to]
     * \param values      The tensor to [read|write]
     */
    void tensor(const std::string& key, ::tensor values) const;

    /*! \brief Return a subset of the current CheckpointData
     *
     * \tparam operation  Whether we are reading or writing
     * \param key         The key to [read|write] the sub data [from|to]
     * \return            A CheckpointData object representing a subset of the current object
     */
    //!{
    CheckpointData subCheckpointData(const std::string& key) const;
    //!}

private:
    //! KV tree read from checkpoint
    const KeyValueTreeObject* inputTree_ = nullptr;

    //! Construct an input checkpoint data object
    explicit CheckpointData(const KeyValueTreeObject& inputTree);

    // Only holders should build
    friend class ReadCheckpointDataHolder;
};

template<>
class CheckpointData<CheckpointDataOperation::Write>
{
public:
    //! \copydoc CheckpointData<CheckpointDataOperation::Read>::scalar
    //! {
    template<typename T>
    std::enable_if_t<IsSerializableType<T>::value, void> scalar(const std::string& key, const T* value);
    template<typename T>
    std::enable_if_t<IsSerializableEnum<T>::value, void> enumScalar(const std::string& key, const T* value);
    //! }

    //! \copydoc CheckpointData<CheckpointDataOperation::Read>::arrayRef
    //! {
    // Write ArrayRef of scalar
    template<typename T>
    std::enable_if_t<IsSerializableType<T>::value, void> arrayRef(const std::string& key,
                                                                  ArrayRef<const T>  values);
    // Write ArrayRef of RVec
    void arrayRef(const std::string& key, ArrayRef<const RVec> values);
    //! }

    //! \copydoc CheckpointData<CheckpointDataOperation::Read>::tensor
    void tensor(const std::string& key, const ::tensor values);

    //! \copydoc CheckpointData<CheckpointDataOperation::Read>::subCheckpointData
    CheckpointData subCheckpointData(const std::string& key);

private:
    //! Builder for the tree to be written to checkpoint
    std::optional<KeyValueTreeObjectBuilder> outputTreeBuilder_ = std::nullopt;

    //! Construct an output checkpoint data object
    explicit CheckpointData(KeyValueTreeObjectBuilder&& outputTreeBuilder);

    // Only holders should build
    friend class WriteCheckpointDataHolder;
};

/*! \brief Read a checkpoint version enum variable
 *
 * This reads the checkpoint version from file. The read version is returned.
 *
 * If the read version is more recent than the code version, this throws an error, since
 * we cannot know what has changed in the meantime. Using newer checkpoint files with
 * old code is not a functionality we can offer. Note, however, that since the checkpoint
 * version is saved by module, older checkpoint files of all simulations that don't use
 * that specific module can still be used.
 *
 * Allowing backwards compatibility of files (i.e., reading an older checkpoint file with
 * a newer version of the code) is in the responsibility of the caller module. They can
 * use the returned file checkpoint version to do that:
 *
 *     const auto fileVersion = checkpointVersion(checkpointData, "version", c_currentVersion);
 *     if (fileVersion >= CheckpointVersion::AddedX)
 *     {
 *         checkpointData->scalar("x", &x_));
 *     }
 *
 * @tparam VersionEnum     The type of the checkpoint version enum
 * @param  checkpointData  A reading checkpoint data object
 * @param  key             The key under which the version is saved - also used for error output
 * @param  programVersion  The checkpoint version of the current code
 * @return                 The checkpoint version read from file
 */
template<typename VersionEnum>
VersionEnum checkpointVersion(const ReadCheckpointData* checkpointData,
                              const std::string&        key,
                              const VersionEnum         programVersion)
{
    VersionEnum fileVersion;
    checkpointData->enumScalar(key, &fileVersion);
    if (fileVersion > programVersion)
    {
        throw FileIOError(
                formatString("The checkpoint file contains a %s that is more recent than the "
                             "current program version and is not backward compatible.",
                             key.c_str()));
    }
    return fileVersion;
}

/*! \brief Write the current code checkpoint version enum variable
 *
 * Write the current program checkpoint version to the checkpoint data object.
 * Returns the written checkpoint version to mirror the signature of the reading version.
 *
 * @tparam VersionEnum     The type of the checkpoint version enum
 * @param  checkpointData  A writing checkpoint data object
 * @param  key             The key under which the version is saved
 * @param  programVersion  The checkpoint version of the current code
 * @return                 The checkpoint version written to file
 */
template<typename VersionEnum>
VersionEnum checkpointVersion(WriteCheckpointData* checkpointData,
                              const std::string&   key,
                              const VersionEnum    programVersion)
{
    checkpointData->enumScalar(key, &programVersion);
    return programVersion;
}

inline ReadCheckpointData::CheckpointData(const KeyValueTreeObject& inputTree) :
    inputTree_(&inputTree)
{
}

inline WriteCheckpointData::CheckpointData(KeyValueTreeObjectBuilder&& outputTreeBuilder) :
    outputTreeBuilder_(outputTreeBuilder)
{
}
/*! \libinternal
 * \brief Holder for read checkpoint data
 *
 * A ReadCheckpointDataHolder object is passed to the checkpoint reading
 * functionality, and then passed into the SimulatorBuilder object. It
 * holds the KV-tree read from file and dispatches CheckpointData objects
 * to the checkpoint clients.
 */
class ReadCheckpointDataHolder
{
public:
    //! Check whether a key exists
    [[nodiscard]] bool keyExists(const std::string& key) const;

    //! Return vector of existing keys
    [[nodiscard]] std::vector<std::string> keys() const;

    //! Deserialize serializer content into the CheckpointData object
    void deserialize(ISerializer* serializer);

    /*! \brief Return a subset of the current CheckpointData
     *
     * \param key         The key to [read|write] the sub data [from|to]
     * \return            A CheckpointData object representing a subset of the current object
     */
    [[nodiscard]] ReadCheckpointData checkpointData(const std::string& key) const;

    //! Write the contents of the Checkpoint to file
    void dump(FILE* out) const;

private:
    //! KV-tree read from checkpoint
    KeyValueTreeObject checkpointTree_;
};

/*! \libinternal
 * \brief Holder for write checkpoint data
 *
 * The WriteCheckpointDataHolder object holds the KV-tree builder and
 * dispatches CheckpointData objects to the checkpoint clients to save
 * their respective data. It is then passed to the checkpoint writing
 * functionality.
 */
class WriteCheckpointDataHolder
{
public:
    //! Serialize the content of the CheckpointData object
    void serialize(ISerializer* serializer);

    /*! \brief Return a subset of the current CheckpointData
     *
     * \param key         The key to [read|write] the sub data [from|to]
     * \return            A CheckpointData object representing a subset of the current object
     */
    [[nodiscard]] WriteCheckpointData checkpointData(const std::string& key);

    /*! \brief
     */
    [[nodiscard]] bool empty() const;

private:
    //! KV-tree builder
    KeyValueTreeBuilder outputTreeBuilder_;
    //! Whether any checkpoint data object has been requested
    bool hasCheckpointDataBeenRequested_ = false;
};

// Function definitions - here to avoid template-related linker problems
// doxygen doesn't like these...
//! \cond
template<typename T>
std::enable_if_t<IsSerializableType<T>::value, void> ReadCheckpointData::scalar(const std::string& key,
                                                                                T* value) const
{
    GMX_RELEASE_ASSERT(inputTree_, "No input checkpoint data available.");
    *value = (*inputTree_)[key].cast<T>();
}

template<typename T>
std::enable_if_t<IsSerializableEnum<T>::value, void> ReadCheckpointData::enumScalar(const std::string& key,
                                                                                    T* value) const
{
    GMX_RELEASE_ASSERT(inputTree_, "No input checkpoint data available.");
    std::underlying_type_t<T> castValue;
    castValue = (*inputTree_)[key].cast<std::underlying_type_t<T>>();
    *value    = static_cast<T>(castValue);
}

template<typename T>
inline std::enable_if_t<IsSerializableType<T>::value, void>
WriteCheckpointData::scalar(const std::string& key, const T* value)
{
    GMX_RELEASE_ASSERT(outputTreeBuilder_, "No output checkpoint data available.");
    outputTreeBuilder_->addValue(key, *value);
}

template<typename T>
inline std::enable_if_t<IsSerializableEnum<T>::value, void>
WriteCheckpointData::enumScalar(const std::string& key, const T* value)
{
    GMX_RELEASE_ASSERT(outputTreeBuilder_, "No output checkpoint data available.");
    auto castValue = static_cast<std::underlying_type_t<T>>(*value);
    outputTreeBuilder_->addValue(key, castValue);
}

template<typename T>
inline std::enable_if_t<IsSerializableType<T>::value, void>
ReadCheckpointData::arrayRef(const std::string& key, ArrayRef<T> values) const
{
    GMX_RELEASE_ASSERT(inputTree_, "No input checkpoint data available.");
    GMX_RELEASE_ASSERT(values.size() >= (*inputTree_)[key].asArray().values().size(),
                       "Read vector does not fit in passed ArrayRef.");
    auto outputIt  = values.begin();
    auto inputIt   = (*inputTree_)[key].asArray().values().begin();
    auto outputEnd = values.end();
    auto inputEnd  = (*inputTree_)[key].asArray().values().end();
    for (; outputIt != outputEnd && inputIt != inputEnd; outputIt++, inputIt++)
    {
        *outputIt = inputIt->cast<T>();
    }
}

template<typename T>
inline std::enable_if_t<IsSerializableType<T>::value, void>
WriteCheckpointData::arrayRef(const std::string& key, ArrayRef<const T> values)
{
    GMX_RELEASE_ASSERT(outputTreeBuilder_, "No output checkpoint data available.");
    auto builder = outputTreeBuilder_->addUniformArray<T>(key);
    for (const auto& value : values)
    {
        builder.addValue(value);
    }
}

inline void ReadCheckpointData::arrayRef(const std::string& key, ArrayRef<RVec> values) const
{
    GMX_RELEASE_ASSERT(values.size() >= (*inputTree_)[key].asArray().values().size(),
                       "Read vector does not fit in passed ArrayRef.");
    auto outputIt  = values.begin();
    auto inputIt   = (*inputTree_)[key].asArray().values().begin();
    auto outputEnd = values.end();
    auto inputEnd  = (*inputTree_)[key].asArray().values().end();
    for (; outputIt != outputEnd && inputIt != inputEnd; outputIt++, inputIt++)
    {
        auto storedRVec = inputIt->asObject()["RVec"].asArray().values();
        *outputIt = { storedRVec[XX].cast<real>(), storedRVec[YY].cast<real>(), storedRVec[ZZ].cast<real>() };
    }
}

inline void WriteCheckpointData::arrayRef(const std::string& key, ArrayRef<const RVec> values)
{
    auto builder = outputTreeBuilder_->addObjectArray(key);
    for (const auto& value : values)
    {
        auto subbuilder = builder.addObject();
        subbuilder.addUniformArray("RVec", { value[XX], value[YY], value[ZZ] });
    }
}

inline void ReadCheckpointData::tensor(const std::string& key, ::tensor values) const
{
    auto array     = (*inputTree_)[key].asArray().values();
    values[XX][XX] = array[0].cast<real>();
    values[XX][YY] = array[1].cast<real>();
    values[XX][ZZ] = array[2].cast<real>();
    values[YY][XX] = array[3].cast<real>();
    values[YY][YY] = array[4].cast<real>();
    values[YY][ZZ] = array[5].cast<real>();
    values[ZZ][XX] = array[6].cast<real>();
    values[ZZ][YY] = array[7].cast<real>();
    values[ZZ][ZZ] = array[8].cast<real>();
}

inline void WriteCheckpointData::tensor(const std::string& key, const ::tensor values)
{
    auto builder = outputTreeBuilder_->addUniformArray<real>(key);
    builder.addValue(values[XX][XX]);
    builder.addValue(values[XX][YY]);
    builder.addValue(values[XX][ZZ]);
    builder.addValue(values[YY][XX]);
    builder.addValue(values[YY][YY]);
    builder.addValue(values[YY][ZZ]);
    builder.addValue(values[ZZ][XX]);
    builder.addValue(values[ZZ][YY]);
    builder.addValue(values[ZZ][ZZ]);
}

inline ReadCheckpointData ReadCheckpointData::subCheckpointData(const std::string& key) const
{
    return CheckpointData((*inputTree_)[key].asObject());
}

inline WriteCheckpointData WriteCheckpointData::subCheckpointData(const std::string& key)
{
    return CheckpointData(outputTreeBuilder_->addObject(key));
}

//! \endcond

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_CHECKPOINTDATA_H
