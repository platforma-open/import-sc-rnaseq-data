import {
  BlockArgs,
  BlockOutputs,
  model,
} from '@platforma-open/milaboratories.import-sc-rnaseq-data.model';
import { awaitStableState, blockTest } from '@platforma-sdk/test';
import { blockSpec as samplesAndDataBlockSpec } from '@platforma-open/milaboratories.samples-and-data';
import { BlockArgs as SamplesAndDataBlockArgs } from '@platforma-open/milaboratories.samples-and-data.model';
import { blockSpec as myBlockSpec } from 'this-block';
import { InferBlockState, uniquePlId, wrapOutputs } from '@platforma-sdk/model';

blockTest(
  'simple project',
  { timeout: 100000 },
  async ({ rawPrj: project, ml, helpers, expect }) => {
    const sndBlockId = await project.addBlock('Samples & Data', samplesAndDataBlockSpec);
    const importBlockId = await project.addBlock('Import scRNA-seq Data', myBlockSpec);

    const sample1Id = uniquePlId();
    const metaColumn1Id = uniquePlId();
    const dataset1Id = uniquePlId();

    const csvHandle = await helpers.getLocalFileHandle('./assets/test_counts.csv');

    project.setBlockArgs(sndBlockId, {
      metadata: [
        {
          id: metaColumn1Id,
          label: 'MetaColumn1',
          global: false,
          valueType: 'Long',
          data: {
            [sample1Id]: 2345
          }
        }
      ],
      sampleIds: [sample1Id],
      sampleLabelColumnLabel: 'Sample Name',
      sampleLabels: { [sample1Id]: 'Sample 1' },
      datasets: [
        {
          id: dataset1Id,
          label: 'Dataset 1',
          content: {
            type: 'Xsv',
            xsvType: 'csv',
            gzipped: false,
            data: {
              [sample1Id]: csvHandle
            }
          }
        }
      ],
      h5adFilesToPreprocess: [],
      seuratFilesToPreprocess: []
    } satisfies SamplesAndDataBlockArgs);
    await project.runBlock(sndBlockId);
    await helpers.awaitBlockDone(sndBlockId, 8000);
    const sndBlockState = project.getBlockState(sndBlockId);
    const importBlockState = project.getBlockState(importBlockId);

    const sdnStableState1 = await helpers.awaitBlockDoneAndGetStableBlockState(sndBlockId, 8000);
    expect(sdnStableState1.outputs).toMatchObject({
      fileImports: { ok: true, value: { [csvHandle]: { done: true } } }
    });

    const importStableState1 = (await awaitStableState(
      importBlockState,
      25000
    )) as InferBlockState<typeof model>;

    expect(importStableState1.outputs).toMatchObject({
      datasetOptions: {
        ok: true,
        value: [
          {
            label: 'Dataset 1'
          }
        ]
      }
    });

    const importStableState1Outputs = wrapOutputs(importStableState1.outputs);

    await project.setBlockArgs(importBlockId, {
      datasetRef: importStableState1Outputs.datasetOptions[0].ref,
    } satisfies BlockArgs);

    await project.runBlock(importBlockId);
    const importStableState2 = (await helpers.awaitBlockDoneAndGetStableBlockState(
      importBlockId,
      100000
    )) as InferBlockState<typeof model>;
    const importOutputs = wrapOutputs<BlockOutputs>(importStableState2.outputs);

    expect(importOutputs.rawCountsPf).toBeDefined();
    expect(importOutputs.cellMetricsPf).toBeDefined();
    expect(importOutputs.resultsSummaryPf).toBeDefined();

    const rawCountsPfHandle = importOutputs.rawCountsPf!;
    const rawCountsPfColumnList = await ml.driverKit.pFrameDriver.listColumns(rawCountsPfHandle);
    expect(rawCountsPfColumnList.length).toBeGreaterThan(0);
  }
);
