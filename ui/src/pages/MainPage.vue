<script setup lang="ts">
import { PlBlockPage, PlBtnGhost, PlDropdownRef, PlMaskIcon24, PlSlideModal } from '@platforma-sdk/ui-vue';
import { useApp } from '../app';
import { plRefsEqual, type PlRef } from '@platforma-sdk/model';

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
      <PlDropdownRef
        v-model="app.model.args.datasetRef"
        :options="app.model.outputs.datasetOptions"
        label="Select dataset"
        clearable
        required
        @update:model-value="setDataset"
      />
    </PlSlideModal>
  </PlBlockPage>
</template>
