//
// Copyright (c) 2008-2017 the Urho3D project.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

#include "../Precompiled.h"

#include "../Core/Context.h"
#include "../Graphics/Camera.h"
#include "../Graphics/DebugRenderer.h"
#include "../Graphics/Geometry.h"
#include "../Graphics/IndexBuffer.h"
#include "../Graphics/Material.h"
#include "../Graphics/OcclusionBuffer.h"
#include "../Graphics/OctreeQuery.h"
#include "../Graphics/TerrainFace.h"
#include "../Graphics/TerrainPatch.h"
#include "../Graphics/Terrain.h"
#include "../Graphics/VertexBuffer.h"
#include "../IO/Log.h"
#include "../Scene/Node.h"

#include "../DebugNew.h"

namespace Urho3D
{

static const float LOD_CONSTANT = 1.0f / 150.0f;

extern const char* GEOMETRY_CATEGORY;


// Create each face detached
TerrainFace::TerrainFace(Context* context) :
    Drawable(context, DRAWABLE_GEOMETRY),
    geometry_(new Geometry(context)),
    maxLodGeometry_(new Geometry(context)),
    occlusionGeometry_(new Geometry(context)),
    vertexBuffer_(new VertexBuffer(context)),
    coordinates_(IntVector2::ZERO),
    m_LodLevel(0),
    m_FaceDirection(QuadFaceInvalid)
{
    // Create new geometry
    geometry_->SetVertexBuffer(0, vertexBuffer_);
    maxLodGeometry_->SetVertexBuffer(0, vertexBuffer_);
    occlusionGeometry_->SetVertexBuffer(0, vertexBuffer_);

}


// Create each face detached
void TerrainFace::SetTerrainFace(unsigned int index, QuadFace facedirection, unsigned maxLod, Terrain * parentTerrain )
{
    m_Index=index;
    m_MaxLod=maxLod;
    m_ParentTerrain=parentTerrain;
    m_FaceDirection = facedirection;

    m_CubeFaceSize=parentTerrain->GetCubeSize();
    if(m_CubeFaceSize==0)
        m_CubeFaceSize=2;
}


TerrainFace::~TerrainFace()
{
}

void TerrainFace::RegisterObject(Context* context)
{
    context->RegisterFactory<TerrainFace>();
}

void TerrainFace::ProcessRayQuery(const RayOctreeQuery& query, PODVector<RayQueryResult>& results)
{

}

void TerrainFace::UpdateBatches(const FrameInfo& frame)
{

}

void TerrainFace::UpdateGeometry(const FrameInfo& frame)
{
    if (vertexBuffer_->IsDataLost())
    {
        if (m_ParentTerrain)
            m_ParentTerrain->BuildTerrain(false,false);
        else
            vertexBuffer_->ClearDataLost();
    }
}

UpdateGeometryType TerrainFace::GetUpdateGeometryType()
{
    // Because there is a latency in starting worker thread updates, and the update of terrain patch LOD should not take
    // much time, always update in the main thread
    return UPDATE_MAIN_THREAD;
}

Geometry* TerrainFace::GetLodGeometry(unsigned batchIndex, unsigned level)
{
    if (!level)
        return maxLodGeometry_;
    else
        return geometry_;
}

unsigned TerrainFace::GetNumOccluderTriangles()
{
    // Check that the material is suitable for occlusion (default material always is)
    Material* mat = batches_[0].material_;
    if (mat && !mat->GetOcclusion())
        return 0;
    else
        return occlusionGeometry_->GetIndexCount() / 3;
}

bool TerrainFace::DrawOcclusion(OcclusionBuffer* buffer)
{
    Material* material = batches_[0].material_;
    if (material)
    {
        if (!material->GetOcclusion())
            return true;
        buffer->SetCullMode(material->GetCullMode());
    }
    else
        buffer->SetCullMode(CULL_CCW);

    const unsigned char* vertexData;
    unsigned vertexSize;
    const unsigned char* indexData;
    unsigned indexSize;
    const PODVector<VertexElement>* elements;

    geometry_->GetRawData(vertexData, vertexSize, indexData, indexSize, elements);
    // Check for valid geometry data
    if (!vertexData || !elements || VertexBuffer::GetElementOffset(*elements, TYPE_VECTOR3, SEM_POSITION) != 0)
    {
        URHO3D_LOGINFO("invalid");
        return false;
    }

    // Draw and check for running out of triangles
    buffer->AddTriangles(node_->GetWorldTransform(), vertexData, vertexSize, geometry_->GetVertexStart(),
                         geometry_->GetVertexCount());

    URHO3D_LOGINFO("draw");
}

void TerrainFace::DrawDebugGeometry(DebugRenderer* debug, bool depthTest)
{

}

void TerrainFace::SetOwner(Terrain* terrain)
{

}

void TerrainFace::SetNeighbors(TerrainFace* north, TerrainFace* south, TerrainFace* west, TerrainFace* east)
{

}

void TerrainFace::SetMaterial(Material* material)
{

}

void TerrainFace::SetBoundingBox(const BoundingBox& box)
{

}

void TerrainFace::SetCoordinates(const IntVector2& coordinates)
{
    coordinates_ = coordinates;
}

void TerrainFace::ResetLod()
{
}

Geometry* TerrainFace::GetGeometry() const
{

}

Geometry* TerrainFace::GetMaxLodGeometry() const
{
    return maxLodGeometry_;
}

Geometry* TerrainFace::GetOcclusionGeometry() const
{
    return occlusionGeometry_;
}

VertexBuffer* TerrainFace::GetVertexBuffer() const
{
    return vertexBuffer_;
}

Terrain* TerrainFace::GetOwner() const
{
    return owner_;
}

void TerrainFace::OnWorldBoundingBoxUpdate()
{
}

unsigned TerrainFace::GetCorrectedLodLevel(unsigned lodLevel)
{
    return lodLevel;
}

void TerrainFace::Build()
{
    // if Parent Terrain is set
    if(m_ParentTerrain)
        m_ParentTerrain->BuildFace(this);

}

void TerrainFace::CreateTerrainPatch(Vector3 coordinates, int facesize, unsigned int GridSize)
{
    // if no parent
    if(!m_ParentTerrain)
        return;

    //create a new patch
    TerrainPatch * CreatedTerrainPatch = new TerrainPatch(context_);

    // Created terrain patch
    if(CreatedTerrainPatch)
    {
        // Use three 3d coordinates
        CreatedTerrainPatch->SetCoordinates(coordinates);

        // Set Face Size Should be enougth to recreate
        CreatedTerrainPatch ->SetFaceSize(facesize); // Set Cube Size/1

        CreatedTerrainPatch->SetFaceDirection(m_FaceDirection);

        CreatedTerrainPatch->SetBasePatchSize(GridSize);

        CreatedTerrainPatch->SetTerrain(m_ParentTerrain);

        // Create a new patch
        m_TerrainPatches.Push(CreatedTerrainPatch);


        m_ParentTerrain->BuildPatch(CreatedTerrainPatch);
    }



}
}
