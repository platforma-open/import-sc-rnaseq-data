import type {
  BlockArgs,
  BlockOutputs,
  model,
} from '@platforma-open/milaboratories.import-sc-rnaseq-data.model';
import { awaitStableState, blockTest } from '@platforma-sdk/test';
import { blockSpec as samplesAndDataBlockSpec } from '@platforma-open/milaboratories.samples-and-data';
import type { BlockArgs as SamplesAndDataBlockArgs } from '@platforma-open/milaboratories.samples-and-data.model';
import { blockSpec as myBlockSpec } from 'this-block';
import type { InferBlockState, PlDataTableModel } from '@platforma-sdk/model';
import { uniquePlId, wrapOutputs } from '@platforma-sdk/model';

// eslint-disable-next-line @typescript-eslint/no-explicit-any
async function validateSampleSummary(
  resultsSummaryModel: PlDataTableModel | undefined,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  ml: any,
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  expect: any,
  expectedCellNumber: number,
): Promise<void> {
  expect(resultsSummaryModel).toBeDefined();
  if (!resultsSummaryModel) throw new Error('resultsSummaryPf is undefined');
  expect(resultsSummaryModel.fullTableHandle).toBeDefined();

  const summarySpecs = await ml.driverKit.pFrameDriver.getSpec(resultsSummaryModel.fullTableHandle);
  expect(summarySpecs.length).toBeGreaterThan(0);

  const columnIndices = summarySpecs.map((_: unknown, i: number) => i);
  const summaryData = await ml.driverKit.pFrameDriver.getData(
    resultsSummaryModel.fullTableHandle,
    columnIndices,
  );
  expect(summaryData.length).toBeGreaterThan(0);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const hasData = summaryData.some((col: any) => col.data.length > 0);
  expect(hasData).toBe(true);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const nCellsColumnIndex = summarySpecs.findIndex(
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (spec: any) => spec.type === 'column' && spec.spec.name === 'pl7.app/rna-seq/nCells',
  );
  expect(nCellsColumnIndex).toBeGreaterThanOrEqual(0);
  const nCellsData = summaryData[nCellsColumnIndex];
  expect(nCellsData).toBeDefined();
  expect(nCellsData.data.length).toBeGreaterThan(0);
  expect(nCellsData.data[0]).toBe(expectedCellNumber);
}

blockTest(
  'counts csv import',
  { timeout: 100000 },
  async ({ rawPrj: project, ml, helpers, expect }) => {
    const sndBlockId = await project.addBlock('Samples & Data', samplesAndDataBlockSpec);
    const importBlockId = await project.addBlock('Import scRNA-seq Data', myBlockSpec);

    const sample1Id = uniquePlId();
    const metaColumn1Id = uniquePlId();
    const dataset1Id = uniquePlId();

    const csvHandle = await helpers.getLocalFileHandle('./assets/test_counts.csv');

    await project.setBlockArgs(sndBlockId, {
      metadata: [
        {
          id: metaColumn1Id,
          label: 'MetaColumn1',
          global: false,
          valueType: 'Long',
          data: {
            [sample1Id]: 2345,
          },
        },
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
              [sample1Id]: csvHandle,
            },
          },
        },
      ],
    } satisfies SamplesAndDataBlockArgs);
    await project.runBlock(sndBlockId);
    await helpers.awaitBlockDone(sndBlockId, 8000);
    const importBlockState = project.getBlockState(importBlockId);

    const sdnStableState1 = await helpers.awaitBlockDoneAndGetStableBlockState(sndBlockId, 8000);
    expect(sdnStableState1.outputs).toMatchObject({
      fileImports: { ok: true, value: { [csvHandle]: { done: true } } },
    });

    const importStableState1 = (await awaitStableState(
      importBlockState,
      25000,
    )) as InferBlockState<typeof model>;

    expect(importStableState1.outputs).toMatchObject({
      datasetOptions: {
        ok: true,
        value: [
          {
            label: 'Dataset 1',
          },
        ],
      },
    });

    const importStableState1Outputs = wrapOutputs(importStableState1.outputs);

    await project.setBlockArgs(importBlockId, {
      datasetRef: importStableState1Outputs.datasetOptions[0].ref,
    } satisfies BlockArgs);

    await project.runBlock(importBlockId);
    const importStableState2 = (await helpers.awaitBlockDoneAndGetStableBlockState(
      importBlockId, 100000));

    const importOutputs = wrapOutputs(importStableState2.outputs as BlockOutputs);

    expect(importOutputs.cellMetricsPf).toBeDefined();
    expect(importOutputs.resultsSummaryPf).toBeDefined();

    await validateSampleSummary(importOutputs.resultsSummaryPf, ml, expect, 5);
  },
);

blockTest(
  'mtx import',
  { timeout: 100000 },
  async ({ rawPrj: project, ml, helpers, expect }) => {
    const sndBlockId = await project.addBlock('Samples & Data', samplesAndDataBlockSpec);
    const importBlockId = await project.addBlock('Import scRNA-seq Data', myBlockSpec);

    const sample1Id = uniquePlId();
    const metaColumn1Id = uniquePlId();
    const dataset1Id = uniquePlId();

    const matrixHandle = await helpers.getLocalFileHandle('./assets/test_matrix.mtx.gz');
    const featuresHandle = await helpers.getLocalFileHandle('./assets/test_features.tsv.gz');
    const barcodesHandle = await helpers.getLocalFileHandle('./assets/test_barcodes.tsv.gz');

    await project.setBlockArgs(sndBlockId, {
      metadata: [
        {
          id: metaColumn1Id,
          label: 'MetaColumn1',
          global: false,
          valueType: 'Long',
          data: {
            [sample1Id]: 2345,
          },
        },
      ],
      sampleIds: [sample1Id],
      sampleLabelColumnLabel: 'Sample Name',
      sampleLabels: { [sample1Id]: 'Sample 1' },
      datasets: [
        {
          id: dataset1Id,
          label: 'Dataset 1',
          content: {
            type: 'CellRangerMTX',
            gzipped: true,
            data: {
              [sample1Id]: {
                'matrix.mtx': matrixHandle,
                'features.tsv': featuresHandle,
                'barcodes.tsv': barcodesHandle,
              },
            },
          },
        },
      ] as unknown as SamplesAndDataBlockArgs['datasets'],
    } satisfies SamplesAndDataBlockArgs);
    await project.runBlock(sndBlockId);
    await helpers.awaitBlockDone(sndBlockId, 8000);
    const importBlockState = project.getBlockState(importBlockId);

    const sdnStableState1 = await helpers.awaitBlockDoneAndGetStableBlockState(sndBlockId, 8000);
    expect(sdnStableState1.outputs).toMatchObject({
      fileImports: {
        ok: true,
        value: {
          [matrixHandle]: { done: true },
          [featuresHandle]: { done: true },
          [barcodesHandle]: { done: true },
        },
      },
    });

    const importStableState1 = (await awaitStableState(
      importBlockState,
      25000,
    )) as InferBlockState<typeof model>;

    expect(importStableState1.outputs).toMatchObject({
      datasetOptions: {
        ok: true,
        value: [
          {
            label: 'Dataset 1',
          },
        ],
      },
    });

    const importStableState1Outputs = wrapOutputs(importStableState1.outputs);

    await project.setBlockArgs(importBlockId, {
      datasetRef: importStableState1Outputs.datasetOptions[0].ref,
    } satisfies BlockArgs);

    await project.runBlock(importBlockId);
    const importStableState2 = (await helpers.awaitBlockDoneAndGetStableBlockState(
      importBlockId, 100000));

    const importOutputs = wrapOutputs(importStableState2.outputs as BlockOutputs);

    expect(importOutputs.cellMetricsPf).toBeDefined();
    expect(importOutputs.resultsSummaryPf).toBeDefined();

    await validateSampleSummary(importOutputs.resultsSummaryPf, ml, expect, 5);
  },
);
