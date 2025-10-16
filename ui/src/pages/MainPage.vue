<script setup lang="ts">
import { PFrameImpl, plRefsEqual, type PlRef } from '@platforma-sdk/model';
import {
  PlAgDataTableV2,
  PlAlert,
  PlBlockPage,
  PlBtnGhost,
  PlBtnGroup,
  PlDropdownRef,
  PlMaskIcon24,
  PlSlideModal,
  usePlDataTableSettingsV2,
  useWatchFetch,
} from '@platforma-sdk/ui-vue';
import { useApp } from '../app';

const app = useApp();

function onModalUpdate(val: boolean) {
  const mustStayOpen = app.model.args.datasetRef === undefined;
  if (mustStayOpen) {
    app.model.ui.settingsOpen = true;
    return;
  }
  app.model.ui.settingsOpen = val;
}

const setDataset = (datasetRef: PlRef | undefined) => {
  app.model.args.datasetRef = datasetRef;
  if (datasetRef)
    app.model.ui.title = 'Import GEX Data - ' + app.model.outputs.datasetOptions?.find((o) => plRefsEqual(o.ref, datasetRef))?.label;
};

const setMatrixDataset = (matrixFileRef: PlRef | undefined) => {
  app.model.args.matrixFileRef = matrixFileRef;
  if (matrixFileRef)
    app.model.ui.title = 'Import GEX Data - ' + app.model.outputs.datasetOptions?.find((o) => plRefsEqual(o.ref, matrixFileRef))?.label;
};

const setBarcodesDataset = (barcodesFileRef: PlRef | undefined) => {
  app.model.args.barcodesFileRef = barcodesFileRef;
};

const setGenesDataset = (genesFileRef: PlRef | undefined) => {
  app.model.args.genesFileRef = genesFileRef;
};

const tableSettings = usePlDataTableSettingsV2({
  model: () => app.model.outputs.resultsSummaryPf,
});

// Get error logs
const errorLogs = useWatchFetch(() => app.model.outputs.errorLog, async (pframeHandle) => {
  if (!pframeHandle) {
    // don't allow user to hit run right after changing dataset, wait for check-format to finish
    app.model.ui.allowRun = false;
    return undefined;
  }
  // Get ID of first pcolumn in the pframe (the only one we will access)
  const pFrame = new PFrameImpl(pframeHandle);
  const list = await pFrame.listColumns();
  const id = list?.[0].columnId;
  if (!id) {
    app.model.ui.allowRun = true;
    return undefined;
  }
  // Get unique values of that first pcolumn
  const response = await pFrame.getUniqueValues({ columnId: id, filters: [], limit: 1000000 });
  if (!response) {
    app.model.ui.allowRun = true;
    return undefined;
  }
  if (response.values.data.length === 0) {
    app.model.ui.allowRun = true;
    return undefined;
  }
  app.model.ui.allowRun = false;
  return response.values.data.join('\n');
});

</script>

<template>
  <PlBlockPage>
    <template #title>
      {{ app.model.ui.title }}
    </template>
    <template #append>
      <PlBtnGhost @click.stop="() => (app.model.ui.settingsOpen = true)">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>

    <PlSlideModal :model-value="app.model.ui.settingsOpen" @update:model-value="onModalUpdate">
      <template #title>Settings</template>

      <PlBtnGroup
        v-model="app.model.args.importMode"
        :options="[
          { value: 'csv', label: 'From CSV' },
          { value: 'mtx', label: 'From MTX' },
        ]"
        label="Importing format"
        tooltip="Select the input format for the gene expression data. Can be either a CSV/TSV file or matrix, barcodes and genes files."
      />

      <template v-if="app.model.args.importMode === 'csv'">
        <PlDropdownRef
          v-model="app.model.args.datasetRef"
          :options="app.model.outputs.datasetOptions"
          label="Select dataset"
          clearable
          required
          @update:model-value="setDataset"
        />
      </template>

      <template v-if="app.model.args.importMode === 'mtx'">
        <PlDropdownRef
          v-model="app.model.args.matrixFileRef"
          :options="app.model.outputs.matrixFileOptions"
          label="Select dataset with matrix files"
          clearable
          required
          @update:model-value="setMatrixDataset"
        />
        <PlDropdownRef
          v-model="app.model.args.barcodesFileRef"
          :options="app.model.outputs.barcodesFileOptions"
          label="Select dataset with barcodes files"
          clearable
          required
          @update:model-value="setBarcodesDataset"
        />
        <PlDropdownRef
          v-model="app.model.args.genesFileRef"
          :options="app.model.outputs.barcodesFileOptions"
          label="Select dataset with genes files"
          clearable
          required
          @update:model-value="setGenesDataset"
        />
      </template>

      <PlAlert v-if="errorLogs.value !== undefined" type="warn" icon>
        {{ errorLogs.value }}
      </PlAlert>
      <PlAlert v-if="app.model.outputs.runningPrerun" type="info" icon>
        "Checking input files' format..."
      </PlAlert>
    </PlSlideModal>
    <PlAgDataTableV2 v-model="app.model.ui.tableState" :settings="tableSettings" show-export-button />
  </PlBlockPage>
</template>

<style scoped>
</style>
